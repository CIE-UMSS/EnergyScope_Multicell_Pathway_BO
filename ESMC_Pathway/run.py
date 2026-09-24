#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ESMC Pathway - main run script
==============================

Runs the ESMC Pathway model (EnergyScope Multi-Cell, 2015-2050) with a
rolling horizon (myopic) approach:
- The pathway is split into windows; each window's decisions are fixed for the next.
- A single window covering all years gives perfect foresight.
- Results are saved to _Results.pkl; use --skip-optimize to only redo post-processing.

Usage
-----
# Full run (optimization + post-processing):
python run.py

# Custom window configuration:
python run.py --window-size 15 --overlap 5

# Skip optimization entirely — load _Results.pkl and regenerate
# all CSV outputs and plots without re-solving:
python run.py --skip-optimize

# Force re-optimization even if _Results.pkl already exists:
python run.py --force-optimize

# Skip TD generation (use existing TD files):
python run.py --skip-td-generation

# Skip .dat regeneration from CSV (use existing .dat files):
python run.py --skip-dat-generation

Arguments
---------
--window-size : int, default=N_YEARS_WND
    Years per optimization window (multiple of 5).

--overlap : int, default=N_YEARS_OVERLAP
    Overlap between consecutive windows in years (multiple of 5)

--skip-optimize : flag
    Skip AMPL optimization and load existing _Results.pkl.
    Useful for fast iteration on post-processing and plots.
    Requires a prior full run to have saved _Results.pkl.

--force-optimize : flag
    Force re-optimization even if _Results.pkl already exists.

--skip-td-generation : flag
    Skip typical day generation (use existing TD files).

--skip-dat-generation : flag
    Skip .dat regeneration from CSV. Use only when you are sure
    the .dat files are already up-to-date.

--force-regenerate : flag
    Force regeneration of all data files.

--td-year : str, default='2050'
    Reference year for TD generation.

--td-algo : str, default='read'
    TD algorithm: kmedoid (generate new) or read (use existing).

--n-seg : int, default=N_SEG
    Intra-day segments per typical day (24 = hourly).

Author: Pablo Jimenez Zabalaga
Based on EnergyScope Multi-Cell (P. Thiran et al.) and the EnergyScope
Pathway rolling horizon implementation.
"""

import numpy as np
from pathlib import Path
import sys
import pandas as pd
import logging
import time
import argparse
from datetime import datetime
import traceback

# ============================================================================
# PROJECT PATH CONFIGURATION
# ============================================================================

# Repository root: folder containing the esmc package (works from the root or from scripts/)
PROJECT_PATH = next(p for p in Path(__file__).resolve().parents if (p / 'esmc').is_dir())
sys.path.insert(0, str(PROJECT_PATH))

# ============================================================================
# IMPORTS (after path configuration; importing esmc also sets up logging)
# ============================================================================

from esmc.utils.esmcpt import EsmcPathway
from esmc.common import bo_country_code, CSV_SEPARATOR

from esmc.utils.esmc_preprocessor import EsmcPreProcessor
from esmc.utils.esmc_collector import EsmcCollector
from esmc.preprocessing.generate_year_dat_files import generate_all_year_dat_files

# ============================================================================
# MODEL CONFIGURATION
# ============================================================================

# ======================
# ROLLING HORIZON CONFIG
# ======================
N_YEARS_WND     = 10   # Window size in years (multiple of 5)
N_YEARS_OVERLAP = 5    # Overlap between windows in years (multiple of 5)

# ======================
# CASE STUDY CONFIG
# ======================
PATHWAY_CASE = 'Six_reg_NZE_2015_2050_SEG_24_hour_twelveTD'

YEARS = ['2015', '2021', '2025', '2030', '2035', '2040', '2045', '2050']

PHASES = {
    '2015_2021': {'start': 'YEAR_2015', 'stop': 'YEAR_2021'},
    '2021_2025': {'start': 'YEAR_2021', 'stop': 'YEAR_2025'},
    '2025_2030': {'start': 'YEAR_2025', 'stop': 'YEAR_2030'},
    '2030_2035': {'start': 'YEAR_2030', 'stop': 'YEAR_2035'},
    '2035_2040': {'start': 'YEAR_2035', 'stop': 'YEAR_2040'},
    '2040_2045': {'start': 'YEAR_2040', 'stop': 'YEAR_2045'},
    '2045_2050': {'start': 'YEAR_2045', 'stop': 'YEAR_2050'}
}

PATHWAY_PARAMS = {
    'limit_LT_renovation':         0.33,
    'limit_pass_mob_changes':      0.50,
    'limit_freight_changes':       0.50,
    'limit_HT_renovation':         0.33,
    'limit_cooking_changes':       0.33,
    'limit_mech_comm_changes':     0.33,
    'limit_mech_mov_agr_changes':  0.33,
    'limit_mech_fix_agr_changes':  0.33,
    'limit_mech_min_changes':      0.33,
    'limit_mech_fish_changes':     0.33
}

NBR_TDS          = 12  # number of typical days
AMPL_PATH        = None
GWP_LIMIT_OVERALL = True
RE_SHARE_PRIMARY  = None
F_PERC            = True
# Model: ESMC_model_AMPL_pathway_SEG.mod (storage level chained over the year,
# with optional intra-day segmentation)
N_SEG             = 24  # intra-day segments per TD (24 = hourly, <24 = segmented with tsam)
SAVE_HOURLY = ['Resources', 'Exchanges', 'Assets', 'Storage', 'Curt']


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def print_header(window_size: int, overlap: int, skip_optimize: bool = False):
    """Print script header with configuration"""
    print(f"\n{'='*70}")
    print(f"ESMC PATHWAY - ROLLING HORIZON OPTIMIZATION")
    print(f"{'='*70}")
    print(f"Configuration:")
    print(f"  Window size : {window_size} years")
    print(f"  Overlap     : {overlap} years")
    print(f"  Typical days: {NBR_TDS}")
    print(f"  Years       : {YEARS[0]} to {YEARS[-1]}")
    if skip_optimize:
        print(f"\n  [FAST MODE] Skipping optimization — loading from _Results.pkl")
    else:
        print(f"\n  [FULL MODE] Running optimization + post-processing")
    print(f"{'='*70}\n")


def print_window_header(window_idx: int, n_windows: int,
                        years: list, phases: list):
    """Print header for each optimization window"""
    print(f"\n{'='*60}")
    print(f"WINDOW {window_idx + 1} of {n_windows}")
    print(f"{'='*60}")
    print(f"Years : {years}")
    print(f"Phases: {phases}")
    print(f"{'-'*60}")


def print_final_summary(total_time: float, solve_time: float,
                        n_windows: int, output_dir: Path,
                        skip_optimize: bool = False):
    """Print final execution summary"""
    print(f"\n{'='*70}")
    if skip_optimize:
        print(f"ROLLING HORIZON — POST-PROCESSING COMPLETED (Fast Mode)")
        print(f"[FAST] Loaded _Results.pkl — no re-solving needed!")
    else:
        print(f"ROLLING HORIZON OPTIMIZATION COMPLETED")
        print(f"Summary:")
        print(f"  Windows optimized     : {n_windows}")
        print(f"  Total solve time      : {solve_time/60:.1f} minutes")
        if n_windows > 0:
            print(f"  Avg time / window     : {solve_time/60/n_windows:.1f} minutes")
    print(f"{'='*70}")
    print(f"Total execution time : {total_time/60:.1f} minutes")
    print(f"Results directory    : {output_dir}")
    if not skip_optimize:
        print(f"\n[TIP] Next time, skip re-solving with:")
        print(f"   python {sys.argv[0]} --skip-optimize")
    print(f"\nOutputs:")
    print(f"  Collector results  : outputs/_Results.pkl  +  outputs/*.csv")
    print(f"  New/Old/Decom      : outputs/New_old_decom_graphs/")
    print(f"  Investment plots   : outputs/Cost_Inv_Phase/")
    print(f"  Operation plots    : outputs/Cost_Op_Phase/")
    print(f"\n{'='*70}")


# ============================================================================
# MAIN FUNCTION
# ============================================================================

def main():
    """Main rolling horizon optimization routine"""

    # ========================================================================
    # ARGUMENT PARSING
    # ========================================================================

    parser = argparse.ArgumentParser(
        description='Rolling Horizon Optimization for ESMC Pathway Model',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run.py                          # full run
  python run.py --window-size 15         # custom window
  python run.py --skip-optimize          # skip re-solving
  python run.py --force-optimize         # force re-solve
  python run.py --skip-td-generation     # reuse TD files
        """
    )

    parser.add_argument('--window-size',   type=int,  default=N_YEARS_WND,
                        help=f'Years per optimization window (default: {N_YEARS_WND})')
    parser.add_argument('--overlap',       type=int,  default=N_YEARS_OVERLAP,
                        help=f'Overlap between windows (default: {N_YEARS_OVERLAP})')
    parser.add_argument('--skip-optimize', action='store_true',
                        help='Skip AMPL optimization; load results from _Results.pkl')
    parser.add_argument('--force-optimize', action='store_true',
                        help='Force re-optimization even if _Results.pkl exists')
    parser.add_argument('--skip-td-generation', action='store_true',
                        help='Skip TD file generation (use existing TD files)')
    parser.add_argument('--skip-dat-generation', action='store_true',
                        help='Skip .dat regeneration from CSV')
    parser.add_argument('--force-regenerate', action='store_true',
                        help='Force regeneration of all data files')
    parser.add_argument('--td-year',  type=str, default='2050',
                        help='Reference year for TD generation (default: 2050)')
    parser.add_argument('--td-algo',  type=str, default='read',
                        choices=['kmedoid', 'read'],
                        help='TD algorithm (default: read)')
    parser.add_argument('--n-seg', type=int, default=N_SEG,
                        help=f'Intra-day segments per TD (default: {N_SEG}; 24 = hourly)')

    args = parser.parse_args()

    if args.skip_optimize and args.force_optimize:
        print("[ERROR] Cannot use both --skip-optimize and --force-optimize")
        return

    n_year_opti    = args.window_size
    n_year_overlap = args.overlap
    skip_optimize  = args.skip_optimize

    print_header(n_year_opti, n_year_overlap, skip_optimize=skip_optimize)
    start_time = time.time()

    # ========================================================================
    # STEP 1: CREATE PATHWAY MODEL
    # ========================================================================

    print("=== Step 1: Creating Pathway Model ===")

    config = {
        'regions_names':    bo_country_code,
        'case_study':       f'{PATHWAY_CASE}_{n_year_opti}y_{n_year_overlap}y',
        'gwp_limit_overall': GWP_LIMIT_OVERALL,
        're_share_primary':  RE_SHARE_PRIMARY,
        'f_perc':            F_PERC,
        'space_id':          '_'.join(bo_country_code)
    }

    years_config = {
        'years':           YEARS,
        'phases':          PHASES,
        'pathway_params':  PATHWAY_PARAMS
    }

    try:
        pathway_model = EsmcPathway(config, years_config, nbr_td=NBR_TDS)
        pathway_model.n_seg = args.n_seg
        print(f'  [MODEL] Intra-day segments: N_SEG={args.n_seg}')
        print(f"[OK] Pathway model created: {config['case_study']}")
        print(f"  Regions          : {len(config['regions_names'])}")
        print(f"  Output directory : {pathway_model.cs_dir}")
    except Exception as e:
        print(f"[ERROR] Failed to create pathway model: {e}")
        traceback.print_exc()
        return

    # Canonical path for the collector pickle
    results_pkl = pathway_model.cs_dir / 'outputs' / '_Results.pkl'

    # Info about existing pkl
    if results_pkl.exists():
        size_mb = results_pkl.stat().st_size / 1e6
        mtime   = datetime.fromtimestamp(results_pkl.stat().st_mtime)
        print(f"\n[INFO] Existing _Results.pkl found:")
        print(f"  Size    : {size_mb:.1f} MB")
        print(f"  Modified: {mtime.strftime('%Y-%m-%d %H:%M')}")
        if not skip_optimize and not args.force_optimize:
            print(f"  [TIP] Use --skip-optimize to skip re-solving")
    else:
        if skip_optimize:
            print(f"\n[ERROR] --skip-optimize requested but no _Results.pkl found at:")
            print(f"  {results_pkl}")
            print(f"  Run a full optimization first.")
            return

    # ========================================================================
    # STEP 2: REGENERATE .DAT FILES FROM CSV DATA
    # Runs before TD generation (also with --skip-optimize).
    # ========================================================================

    if not args.skip_dat_generation:
        print("\n=== Step 2: Regenerating .dat files from CSV data ===")
        print("[INFO] Compiling year .dat files from source CSVs in Data/...")

        try:
            dat_success = generate_all_year_dat_files(
                data_dir=pathway_model.data_dir,
                regions_names=list(config['regions_names']),
                years=YEARS,
                verbose=True
            )
        except Exception as e:
            print(f"[ERROR] Failed to regenerate .dat files: {e}")
            traceback.print_exc()
            dat_success = False

        if not dat_success:
            print("[ERROR] .dat file generation failed or incomplete.")
            print("  Use --skip-dat-generation to proceed with existing files.")
            return
        else:
            print("[OK] All year .dat files regenerated from CSV data\n")
    else:
        print("\n=== Step 2: Skipping .dat regeneration (--skip-dat-generation) ===\n")

    # ========================================================================
    # STEP 3: GENERATE TYPICAL DAYS (if needed)
    # ========================================================================

    if not args.skip_td_generation:
        print("\n=== Step 3: Generating Typical Days ===")
        try:
            pathway_model.read_data_indep(year=args.td_year)
            pathway_model.init_regions(year=args.td_year)
            td_algo = 'kmedoid' if args.force_regenerate else args.td_algo
            pathway_model.init_ta(algo=td_algo, ampl_path=AMPL_PATH)
            pathway_model.print_td_data()
            print(f"[OK] Typical days generated ({NBR_TDS} TDs)")
        except Exception as e:
            print(f"[ERROR] Failed to generate typical days: {e}")
            traceback.print_exc()
            return
    else:
        print("\n=== Step 3: Skipping TD Generation (using existing files) ===")

    # ========================================================================
    # STEP 4: PROCESS PATHWAY DATA FILES
    # ========================================================================

    print("\n=== Step 4: Processing Pathway Data ===")
    try:
        pathway_model.process_pathway_data()
        print("[OK] Pathway data files processed")
    except Exception as e:
        print(f"[ERROR] Failed to process pathway data: {e}")
        traceback.print_exc()
        return

    # ========================================================================
    # STEP 5: ROLLING HORIZON LOOP  –OR–  LOAD FROM _Results.pkl
    # ========================================================================

    total_solve_time = 0
    window_times     = []
    n_windows        = 0
    esmc_collector   = None

    if not skip_optimize:

        # ====================================================================
        # STEP 5a: INITIALIZE ROLLING HORIZON
        # ====================================================================

        print("\n=== Step 5: Initializing Rolling Horizon ===")
        try:
            pathway_model.clean_history()

            esmc_pre = EsmcPreProcessor(
                pathway_model,
                n_years_wnd=n_year_opti,
                n_years_overlap=n_year_overlap
            )
            esmc_pre.print_window_summary()

            output_file = results_pkl
            expl_text   = (f'Rolling horizon: {n_year_opti}y window, '
                           f'{n_year_overlap}y overlap, {NBR_TDS} TDs')

            esmc_collector = EsmcCollector(esmc_pre, output_file, expl_text)
            n_windows      = esmc_pre.get_num_windows()

            print(f"\n[OK] Rolling horizon initialized")
            print(f"  Number of windows: {n_windows}")
            print(f"  Results will be saved to: {output_file}")

        except Exception as e:
            print(f"[ERROR] Failed to initialize rolling horizon: {e}")
            traceback.print_exc()
            return

        # ====================================================================
        # STEP 5b: ROLLING HORIZON OPTIMIZATION LOOP
        # ====================================================================

        print("\n=== Step 6: Rolling Horizon Optimization Loop ===")

        for i in range(n_windows):
            window_start = time.time()

            info = esmc_pre.get_window_info(i)
            print_window_header(i, n_windows, info['years'], info['phases'])

            try:
                # Write sequential optimization file
                print(f"  Writing PES_seq_opti.dat...")
                curr_years_wnd = esmc_pre.write_seq_opti(i)

                # Update remaining years file
                print(f"  Updating remaining years...")
                esmc_pre.update_remaining_years(i)

                # Set up AMPL model (window index selects CPLEX options)
                print(f"  Setting up AMPL model...")
                pathway_model._window_idx = i
                if i == 0:
                    pathway_model.set_esom_pathway(ampl_path=AMPL_PATH, include_fix=False)
                else:
                    pathway_model.reload_with_fix(ampl_path=AMPL_PATH)

                # Solve optimization
                print(f"  Solving optimization...")
                solve_start = time.time()
                pathway_model.solve_esom_pathway()
                solve_time   = time.time() - solve_start
                total_solve_time += solve_time
                print(f"  [OK] Solved in {solve_time/60:.1f} minutes")

                try:
                    solve_result = pathway_model.esom.ampl.getValue("solve_result")
                    print(f"  Solver status: {solve_result}")
                except Exception:
                    pass

                # Extract year-indexed results into pathway_model.results
                print(f"  Extracting year-indexed results...")
                pathway_model.get_total_cost_pathway()
                pathway_model.get_cost_breakdown_pathway()
                pathway_model.get_gwp_breakdown_pathway()
                pathway_model.get_resources_and_exchanges_pathway(save_hourly=SAVE_HOURLY)
                pathway_model.get_assets_pathway(save_hourly=SAVE_HOURLY)
                pathway_model.get_year_balance_pathway()
                pathway_model.get_curt_pathway(save_hourly=SAVE_HOURLY)
                print(f"  [OK] Extracted {len(pathway_model.results)} year-indexed result types")

                # Initialize collector on first window
                if i == 0:
                    esmc_collector.init_storage(pathway_model)

                # Exclude overlap year from year-indexed collection
                if i > 0 and esmc_pre.year_to_rm in curr_years_wnd:
                    curr_years_wnd = [y for y in curr_years_wnd
                                      if y != esmc_pre.year_to_rm]
                    print(f"  Excluding overlap year: {esmc_pre.year_to_rm}")

                # Accumulate year-indexed results in collector
                print(f"  Accumulating year-indexed results in collector...")
                esmc_collector.update_storage(pathway_model, curr_years_wnd, i)

                # Extract phase-indexed variables while AMPL is still alive
                # (New_old_decom, C_inv_phase_tech, C_op_phase_tech, C_op_phase_res)
                print(f"  Extracting phase-indexed results from AMPL...")
                phases_wnd = info['phases']  # phases optimized in this window
                # Window 0: add initial phase '2015_2021' to capture pre-existing capacity
                is_first_window = (i == 0)
                if is_first_window:
                    phases_wnd = ['2015_2021'] + phases_wnd
                esmc_collector.extract_phase_results_from_ampl(
                    pathway_model, phases_wnd,
                    is_first_window=is_first_window
                )
                print(f"  [OK] Phase-indexed results accumulated")

                # Fix solution for next window (except last)
                if i < n_windows - 1:
                    print(f"  Fixing solution for next window...")
                    pathway_model.set_init_sol()

                window_time = time.time() - window_start
                window_times.append(window_time)
                print(f"\n  Window {i+1} completed in {window_time/60:.1f} minutes")

            except Exception as e:
                print(f"[ERROR] Failed on window {i+1}: {e}")
                traceback.print_exc()
                print(f"  Attempting to continue with next window...")
                continue

        # ====================================================================
        # STEP 7: FINALIZE AND SAVE COLLECTOR (_Results.pkl, used by --skip-optimize)
        # ====================================================================

        print("\n=== Step 7: Finalizing and Saving Collector Results ===")
        try:
            esmc_collector.clean_collector()
            esmc_collector.aggregate_regions()

            # Save _Results.pkl  (the single source of truth for skip-optimize)
            esmc_collector.save_results()
            print(f"[OK] _Results.pkl saved to: {results_pkl}")

            # Save CSV files (regional + aggregated)
            csv_dir = pathway_model.cs_dir / 'outputs'
            esmc_collector.save_csv(csv_dir)
            esmc_collector.print_summary()

            print("[OK] Collector results finalized")

        except Exception as e:
            print(f"[ERROR] Failed to finalize collector results: {e}")
            traceback.print_exc()

        # Close last AMPL instance
        try:
            pathway_model.esom.ampl.close()
        except Exception:
            pass

    else:
        # ====================================================================
        # SKIP-OPTIMIZE PATH: LOAD FROM _Results.pkl
        # ====================================================================

        print("\n=== Step 5: Loading Results from _Results.pkl ===")
        try:
            esmc_collector = EsmcCollector.load_from_pkl(results_pkl)
            print(f"[OK] Loaded results from: {results_pkl}")
            esmc_collector.print_summary()
        except Exception as e:
            print(f"[ERROR] Failed to load _Results.pkl: {e}")
            traceback.print_exc()
            return

    # ========================================================================
    # POST-PROCESSING (full run and --skip-optimize)
    # restore_to_pathway_model() loads results and a MockEsom (plots without AMPL)
    # ========================================================================

    print("\n=== Step 8: Restoring Full Pathway Results to Model ===")
    try:
        esmc_collector.restore_to_pathway_model(pathway_model)
        print("[OK] Full pathway results restored to pathway_model")
        print(f"  Available result keys: "
              f"{[k for k, v in pathway_model.results.items() if v is not None]}")
    except Exception as e:
        print(f"[ERROR] Failed to restore results: {e}")
        traceback.print_exc()
        return

    # ========================================================================
    # STEP 9: NEW / OLD / DECOM
    # ========================================================================

    print("\n=== Step 9: Processing New/Old/Decom Results ===")
    output_dir = pathway_model.cs_dir / 'outputs'
    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        new_old_decom_df = pathway_model.results.get('New_old_decom')

        if new_old_decom_df is not None and not new_old_decom_df.empty:
            # Save/re-save the CSV
            csv_path = output_dir / 'New_old_decom.csv'
            new_old_decom_df.to_csv(csv_path, sep=';')
            print(f"  [OK] New_old_decom.csv saved ({len(new_old_decom_df)} rows)")

            # Generate plots
            pathway_model.plot_new_old_decom_by_sector(
                new_old_decom_df=new_old_decom_df,
                output_dir=output_dir
            )
            print("[OK] New/Old/Decom plots saved to outputs/New_old_decom_graphs/")
        else:
            print("[WARN] New_old_decom data not found — skipping plots")

    except Exception as e:
        print(f"[WARN] Error processing New/Old/Decom: {e}")
        traceback.print_exc()

    # ========================================================================
    # STEP 10: COST PLOTS
    # (graph_cost_* use self.esom.get_var() which now reads from MockEsom)
    # ========================================================================

    print("\n=== Step 10: Generating Cost Plots ===")
    try:
        plot_start = time.time()

        print("  Generating investment cost plots...")
        df_inv = pathway_model.graph_cost_inv_phase_tech(plot=True)
        if df_inv is not None:
            print("[OK] Investment plots → outputs/Cost_Inv_Phase/")

        print("  Generating operation cost plots...")
        df_op = pathway_model.graph_cost_op_phase(plot=True)
        if df_op is not None:
            print("[OK] Operation plots → outputs/Cost_Op_Phase/")

        print(f"[OK] Cost plots generated in {time.time() - plot_start:.1f} s")

    except Exception as e:
        print(f"[WARN] Error generating cost plots: {e}")
        traceback.print_exc()

    # ========================================================================
    # FINAL SUMMARY
    # ========================================================================

    total_time = time.time() - start_time
    print_final_summary(
        total_time=total_time,
        solve_time=total_solve_time,
        n_windows=n_windows,
        output_dir=output_dir,
        skip_optimize=skip_optimize
    )

    if not skip_optimize and window_times:
        print(f"\nPer-window timing:")
        for i, t in enumerate(window_times):
            print(f"  Window {i+1}: {t/60:.1f} minutes")
        print(f"  Average : {np.mean(window_times)/60:.1f} minutes")
        print(f"  Std dev : {np.std(window_times)/60:.1f} minutes")

    print(f"\n[COMPLETE] Rolling horizon finished.")


# ============================================================================
# ENTRY POINT
# ============================================================================

if __name__ == "__main__":
    main()