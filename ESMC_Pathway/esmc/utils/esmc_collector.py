#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Results Collector for ESMC Rolling Horizon Optimization
=======================================================

This module handles the collection and aggregation of results across
multiple optimization windows in a rolling horizon framework.

Key responsibilities:
- Initialize storage structures for all result types
- Update storage after each window optimization
- Handle overlap between windows (avoid double-counting)
- Save aggregated results to files

Based on AmplCollector from EnergyScope Pathway (single-region),
adapted for the multi-regional ESMC model.

Author: Pablo Jimenez Zabalaga (adapted for ESMC)
Date: January 2026
"""

from pathlib import Path
import os
import csv
import pickle
import pandas as pd
import numpy as np
from datetime import datetime
from time import time
import logging
from typing import Dict, List, Optional, Any


class EsmcCollector:
    """
    Collector for aggregating results across rolling horizon optimization windows.

    This class manages the collection of results from each optimization window,
    handling the merging of results while avoiding double-counting of overlapping
    periods.

    Parameters
    ----------
    esmc_pre : EsmcPreProcessor
        The preprocessor object with window definitions
    output_file : Path or str
        Path to save the collected results (pickle format)
    expl_text : str
        Description of the case study for documentation

    Attributes
    ----------
    results : dict
        Dictionary of collected results DataFrames
    results_all : dict
        Dictionary of region-aggregated results

    Example
    -------
    >>> from esmc.utils.esmc_preprocessor import EsmcPreProcessor
    >>> from esmc.utils.esmc_collector import EsmcCollector
    >>>
    >>> pre = EsmcPreProcessor(pathway, n_years_wnd=10, n_years_overlap=5)
    >>> collector = EsmcCollector(pre, output_file, "Rolling horizon test")
    >>>
    >>> for i in range(pre.get_num_windows()):
    ...     # ... optimize window ...
    ...     if i == 0:
    ...         collector.init_storage(pathway)
    ...     collector.update_storage(pathway, years, i)
    >>>
    >>> collector.clean_collector()
    >>> collector.save_results()
    """

    def __init__(self, esmc_pre, output_file, expl_text: str):
        self.esmc_pre = esmc_pre
        self.pth_output_all = Path(output_file).parent.parent
        self.output_file = Path(output_file)
        self.expl_text = expl_text

        # Initialize empty result dictionaries
        self.results: Dict[str, Optional[pd.DataFrame]] = {}
        self.results_all: Dict[str, Optional[pd.DataFrame]] = {}
        self.hourly_results: Dict[str, Optional[pd.DataFrame]] = {}

        # Track which windows have been processed
        self._windows_processed: List[int] = []

        logging.info(f'EsmcCollector initialized')
        logging.info(f'  Output file: {self.output_file}')

    def init_storage(self, esmc_pathway) -> None:
        """
        Initialize storage structures based on model configuration.

        This should be called after the first window optimization to set up
        the result storage structures with the correct indices.

        Parameters
        ----------
        esmc_pathway : EsmcPathway
            The pathway model after first window optimization
        """
        # Get full set of years and phases (not just first window)
        Years = ['YEAR_' + y for y in esmc_pathway.years]
        if 'YEAR_2015' in Years:
            Years.remove('YEAR_2015')

        Phases = list(esmc_pathway.phases.keys())
        Regions = esmc_pathway.regions_names

        # Define result categories and their index structures
        # Year-indexed results (primary index is Years)
        year_indexed = [
            'TotalCost', 'TotalGwp', 'Cost_breakdown', 'Gwp_breakdown',
            'Resources', 'Assets', 'Sto_assets', 'Year_balance',
            'Transfer_capacity', 'Exchanges_year', 'Curt',
            'Exch_freight', 'Exch_freight_border'
        ]

        # Phase-indexed results (primary index is Phases)
        phase_indexed = [
            'New_old_decom', 'F_decom', 'C_inv_phase',
            'C_inv_phase_tech', 'C_op_phase_tech', 'C_op_phase_res'
        ]

        # Initialize all result dictionaries
        all_keys = year_indexed + phase_indexed
        self.results = {k: None for k in all_keys}
        self.results_all = {k: None for k in all_keys}

        logging.info(f'Storage initialized for {len(all_keys)} result types')
        logging.info(f'  Year-indexed: {len(year_indexed)}')
        logging.info(f'  Phase-indexed: {len(phase_indexed)}')

    def update_storage(self, esmc_pathway, curr_years_wnd: List[str],
                       window_idx: int) -> None:
        """
        Update storage with results from current optimization window.

        This method extracts results from the pathway model and merges them
        with previously collected results. For overlapping periods, the most
        recent window's results take precedence.

        Parameters
        ----------
        esmc_pathway : EsmcPathway
            The pathway model after window optimization
        curr_years_wnd : list of str
            Years in current window (may be modified to exclude overlap)
        window_idx : int
            Index of current window (0-based)
        """
        # Get phases up to current window (including 2015_2021 initial phase)
        phases_up_to = ['2015_2021'] + self.esmc_pre.phases_up_to[window_idx]
        phases_wnd = self.esmc_pre.phases_opti[window_idx]

        # Track processed window
        self._windows_processed.append(window_idx)

        # Check if model has results
        if not hasattr(esmc_pathway, 'results') or esmc_pathway.results is None:
            logging.warning(f'No results found in pathway model for window {window_idx + 1}')
            return

        n_updated = 0

        for k in self.results:
            # Skip if this result type doesn't exist in model
            if k not in esmc_pathway.results or esmc_pathway.results[k] is None:
                continue

            results_df = esmc_pathway.results[k]

            if results_df is None or (isinstance(results_df, pd.DataFrame) and results_df.empty):
                continue

            # Filter results based on index type
            temp_res = self._filter_results(results_df, k, curr_years_wnd,
                                           phases_up_to, phases_wnd)

            if temp_res is None or temp_res.empty:
                continue

            # Merge with existing results
            if self.results[k] is None:
                self.results[k] = temp_res.copy()
            else:
                # Concatenate and remove duplicates (keep latest)
                combined = pd.concat([self.results[k], temp_res])
                self.results[k] = combined.loc[~combined.index.duplicated(keep='last')]
                self.results[k] = self.results[k].sort_index()

            n_updated += 1

        # ----- Accumulate HOURLY results across windows -----
        # Same year filtering and last-wins dedup as the annual results above.
        n_hourly = 0
        if getattr(esmc_pathway, 'hourly_results', None):
            for k, hdf in esmc_pathway.hourly_results.items():
                if hdf is None or (isinstance(hdf, pd.DataFrame) and hdf.empty):
                    continue
                temp_h = self._filter_results(hdf, k, curr_years_wnd,
                                              phases_up_to, phases_wnd)
                if temp_h is None or temp_h.empty:
                    continue
                if self.hourly_results.get(k) is None:
                    self.hourly_results[k] = temp_h.copy()
                else:
                    combined_h = pd.concat([self.hourly_results[k], temp_h])
                    self.hourly_results[k] = combined_h.loc[
                        ~combined_h.index.duplicated(keep='last')].sort_index()
                n_hourly += 1

        logging.info(f'Updated storage for window {window_idx + 1}: '
                     f'{n_updated} result types, {n_hourly} hourly types')

    def _filter_results(self, results_df: pd.DataFrame, result_key: str,
                       curr_years_wnd: List[str], phases_up_to: List[str],
                       phases_wnd: List[str]) -> Optional[pd.DataFrame]:
        """
        Filter results DataFrame based on result type and current window.

        Parameters
        ----------
        results_df : pd.DataFrame
            Results DataFrame from model
        result_key : str
            Name of the result type
        curr_years_wnd : list
            Years in current window
        phases_up_to : list
            Phases up to current window
        phases_wnd : list
            Phases in current window

        Returns
        -------
        pd.DataFrame or None
            Filtered results, or None if filtering failed
        """
        # Phase-indexed results
        phase_indexed = ['New_old_decom', 'F_decom', 'C_inv_phase',
                        'C_inv_phase_tech', 'C_op_phase_tech', 'C_op_phase_res']

        try:
            if result_key in phase_indexed:
                # Filter by phase
                if 'Phases' in results_df.index.names:
                    return results_df.loc[
                        results_df.index.get_level_values('Phases').isin(phases_up_to), :
                    ]
                elif 'Phase' in results_df.index.names:
                    return results_df.loc[
                        results_df.index.get_level_values('Phase').isin(phases_up_to), :
                    ]
                else:
                    return results_df
            else:
                # Filter by year. curr_years_wnd uses 'YEAR_2021' but some results
                # (e.g. Year_balance) use '2021', so both forms are accepted.
                years_set = set(curr_years_wnd)
                # Add un-prefixed forms: 'YEAR_2021' -> '2021'
                years_set |= {y.replace('YEAR_', '') for y in curr_years_wnd}
                # Add prefixed forms in case any entry is already without prefix
                years_set |= {'YEAR_' + y.replace('YEAR_', '') for y in curr_years_wnd}

                if 'Years' in results_df.index.names:
                    return results_df.loc[
                        results_df.index.get_level_values('Years').isin(years_set), :
                    ]
                elif 'Year' in results_df.index.names:
                    return results_df.loc[
                        results_df.index.get_level_values('Year').isin(years_set), :
                    ]
                else:
                    return results_df

        except Exception as e:
            logging.warning(f'Could not filter {result_key}: {e}')
            return None

    def clean_collector(self) -> None:
        """
        Clean collected results by removing NaN rows.

        This should be called after all windows have been processed.
        """
        n_cleaned = 0

        for k in self.results:
            if self.results[k] is not None and isinstance(self.results[k], pd.DataFrame):
                original_len = len(self.results[k])
                self.results[k].dropna(how='all', inplace=True)
                new_len = len(self.results[k])
                if new_len < original_len:
                    n_cleaned += original_len - new_len

        logging.info(f'Cleaned collector: removed {n_cleaned} NaN rows')

    def aggregate_regions(self) -> None:
        """
        Create region-aggregated versions of results.

        For each result type, creates a version summed across regions
        (where appropriate).
        """
        # Result types that must stay regional-only (never aggregated to a
        # national/root CSV), as in the original single-year ESMC, where
        # Exch_freight_border is stored in self.results but not self.results_all.
        regional_only = {'Exch_freight_border'}

        for k, df in self.results.items():
            if df is None or df.empty:
                self.results_all[k] = None
                continue

            # Skip aggregation for regional-only results: they should only be
            # written under regional_results/, not at the outputs root.
            if k in regional_only:
                self.results_all[k] = None
                continue

            try:
                if 'Regions' in df.index.names:
                    # Sum across regions
                    agg_df = df.groupby(
                        [n for n in df.index.names if n != 'Regions']
                    ).sum()
                    self.results_all[k] = agg_df
                else:
                    # No regions index, use as-is
                    self.results_all[k] = df.copy()
            except Exception as e:
                logging.warning(f'Could not aggregate {k}: {e}')
                self.results_all[k] = None

        logging.info('Created region-aggregated results')

    def save_results(self) -> Path:
        """
        Save collected results to pickle file.

        Also updates the recap CSV file with case study information.

        Returns
        -------
        Path
            Path to the saved pickle file
        """
        case_name = os.path.basename(os.path.normpath(self.output_file.parent))

        # Update recap file
        recap_file = self.pth_output_all / '_Recap.csv'
        t = datetime.now()

        if not recap_file.exists():
            recap_file.parent.mkdir(parents=True, exist_ok=True)
            with open(recap_file, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['Case_study', 'Comment', 'Date_Time', 'N_Windows'])

        # Read existing recap
        try:
            df = pd.read_csv(recap_file)
            if case_name in df.Case_study.values:
                df.loc[df['Case_study'] == case_name, 'Comment'] = self.expl_text
                df.loc[df['Case_study'] == case_name, 'Date_Time'] = t
                df.loc[df['Case_study'] == case_name, 'N_Windows'] = len(self._windows_processed)
                df.to_csv(recap_file, index=False)
            else:
                with open(recap_file, 'a', newline='') as f:
                    writer = csv.writer(f)
                    writer.writerow([case_name, self.expl_text, t, len(self._windows_processed)])
        except Exception as e:
            logging.warning(f'Could not update recap file: {e}')

        # Create output directory if needed
        self.output_file.parent.mkdir(parents=True, exist_ok=True)

        # Prepare data for saving
        save_data = {
            'results': self.results,
            'results_all': self.results_all,
            'metadata': {
                'expl_text': self.expl_text,
                'timestamp': t.isoformat(),
                'windows_processed': self._windows_processed,
                'n_windows': len(self._windows_processed)
            }
        }

        # Save to pickle
        with open(self.output_file, 'wb') as f:
            pickle.dump(save_data, f, protocol=pickle.HIGHEST_PROTOCOL)

        logging.info(f'Results saved to {self.output_file}')
        return self.output_file

    def save_csv(self, output_dir: Path, separator: str = ';') -> None:
        """
        Save all results as CSV files.

        Parameters
        ----------
        output_dir : Path
            Directory to save CSV files
        separator : str, default=';'
            CSV field separator
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Save regional results
        regional_dir = output_dir / 'regional_results'
        regional_dir.mkdir(exist_ok=True)

        n_saved = 0
        for k, df in self.results.items():
            if df is not None and not df.empty:
                output_file = regional_dir / f'{k}.csv'
                df.to_csv(output_file, sep=separator)
                n_saved += 1

        logging.info(f'Saved {n_saved} regional result files to {regional_dir}')

        # Save aggregated results
        if self.results_all:
            n_saved_agg = 0
            for k, df in self.results_all.items():
                if df is not None and not df.empty:
                    output_file = output_dir / f'{k}.csv'
                    df.to_csv(output_file, sep=separator)
                    n_saved_agg += 1

            logging.info(f'Saved {n_saved_agg} aggregated result files to {output_dir}')

        # Save hourly results (written to a dedicated subfolder, if collected)
        if getattr(self, 'hourly_results', None):
            hourly_dir = output_dir / 'hourly_results'
            n_saved_h = 0
            for k, df in self.hourly_results.items():
                if df is not None and not df.empty:
                    hourly_dir.mkdir(exist_ok=True)
                    df.to_csv(hourly_dir / f'{k}.csv', sep=separator)
                    n_saved_h += 1
            if n_saved_h:
                logging.info(f'Saved {n_saved_h} hourly result files to {hourly_dir}')

    def get_summary(self) -> dict:
        """
        Get a summary of collected results.

        Returns
        -------
        dict
            Summary information about collected results
        """
        summary = {
            'windows_processed': len(self._windows_processed),
            'result_types': {},
            'total_rows': 0
        }

        for k, df in self.results.items():
            if df is not None and not df.empty:
                summary['result_types'][k] = {
                    'rows': len(df),
                    'columns': len(df.columns),
                    'index_names': list(df.index.names)
                }
                summary['total_rows'] += len(df)

        return summary


    # =========================================================================
    # PHASE-LEVEL EXTRACTION  (call inside rolling-horizon loop)
    # =========================================================================

    def extract_phase_results_from_ampl(self, pathway_model, phases_wnd: list,
                                         is_first_window: bool = False) -> None:
        """
        Extract phase-indexed variables from AMPL and accumulate in self.results.

        Must be called INSIDE the rolling-horizon loop, after solving each window
        but BEFORE AMPL is closed / reloaded for the next window.

        Extracts and accumulates:
          - New_old_decom  (F_new, F_old, F_decom joined)
          - C_inv_phase_tech
          - C_op_phase_tech
          - C_op_phase_res

        Overlapping phases are de-duplicated using keep='last' (most recent
        window wins), consistent with the year-indexed merge in update_storage.

        Parameters
        ----------
        pathway_model : EsmcPathway
            The pathway model after window optimization (AMPL still live).
        phases_wnd : list of str
            Phases to extract. For window 0 this should include '2015_2021'.
        is_first_window : bool, default False
            If True, also compute phase costs for '2015_2021' manually,
            since AMPL does not compute C_inv_phase_tech / C_op_phase_tech /
            C_op_phase_res for '2015_2021' in rolling horizon mode (it is
            not in PHASE_WND).
        """
        import numpy as np

        ampl = pathway_model.esom.ampl

        # ------------------------------------------------------------------
        # 1. New / Old / Decom
        # ------------------------------------------------------------------
        try:
            # --- F_new ---
            f_new_raw = ampl.getVariable('F_new').getValues().toPandas()
            val_col   = [c for c in f_new_raw.columns if '.val' in str(c)][0]
            f_new_df  = f_new_raw[[val_col]].rename(columns={val_col: 'F_new'})

            n_lvl = f_new_df.index.nlevels
            is_multi = n_lvl == 3
            idx_names = (['Phases', 'Regions', 'Technologies'] if is_multi
                         else ['Phases', 'Technologies'])
            f_new_df.index.names = idx_names

            # --- F_old ---
            f_old_raw = ampl.getVariable('F_old').getValues().toPandas()
            val_col   = [c for c in f_old_raw.columns if '.val' in str(c)][0]
            f_old_df  = f_old_raw[[val_col]].rename(columns={val_col: 'F_old'})
            f_old_df.index.names = idx_names

            # --- F_decom (aggregate over phase_built dimension) ---
            f_decom_raw = ampl.getVariable('F_decom').getValues().toPandas()
            val_col     = [c for c in f_decom_raw.columns if '.val' in str(c)][0]
            f_decom_df  = f_decom_raw[[val_col]].rename(columns={val_col: 'F_decom'})

            if is_multi:
                if f_decom_df.index.nlevels == 4:
                    f_decom_df = f_decom_df.groupby(level=[0, 2, 3]).sum()
                f_decom_df.index.names = ['Phases', 'Regions', 'Technologies']
            else:
                if f_decom_df.index.nlevels == 3:
                    f_decom_df = f_decom_df.groupby(level=[0, 2]).sum()
                f_decom_df.index.names = ['Phases', 'Technologies']

            # --- filter storage technologies ---
            try:
                sto_tech = list(ampl.getSet('STORAGE_TECH').getValues())
                f_decom_df = f_decom_df.loc[
                    ~f_decom_df.index.get_level_values('Technologies').isin(sto_tech), :
                ]
                f_new_df = f_new_df.loc[
                    ~f_new_df.index.get_level_values('Technologies').isin(sto_tech), :
                ]
                f_old_df = f_old_df.loc[
                    ~f_old_df.index.get_level_values('Technologies').isin(sto_tech), :
                ]
            except Exception:
                pass

            # --- join & threshold ---
            nod = f_new_df.join(f_old_df, how='outer').join(f_decom_df, how='outer')
            nod.fillna(0, inplace=True)
            threshold = 1e-5
            nod = nod.mask(nod.abs() < threshold, other=np.nan)

            # --- filter to phases of THIS window only ---
            if 'Phases' in nod.index.names:
                nod = nod.loc[
                    nod.index.get_level_values('Phases').isin(phases_wnd), :
                ]

            # --- accumulate (dedup by phase, keep latest) ---
            prev = self.results.get('New_old_decom')
            if prev is None or (isinstance(prev, pd.DataFrame) and prev.empty):
                self.results['New_old_decom'] = nod
            else:
                combined = pd.concat([prev, nod])
                self.results['New_old_decom'] = combined.loc[
                    ~combined.index.duplicated(keep='last')
                ].sort_index()

            logging.info(f'  [collector] New_old_decom accumulated: '
                         f'{len(self.results["New_old_decom"])} rows')

        except Exception as e:
            logging.warning(f'  [collector] Could not extract New_old_decom: {e}')

        # ------------------------------------------------------------------
        # 2. Phase-cost variables  (C_inv_phase_tech, C_op_phase_tech, C_op_phase_res)
        # ------------------------------------------------------------------
        phase_vars = {
            'C_inv_phase_tech': 'Phase',
            'C_op_phase_tech':  'Phase',
            'C_op_phase_res':   'Phase',
        }

        for var_name, phase_idx_name in phase_vars.items():
            try:
                df = pathway_model.esom.get_var(var_name).copy()
                if df.empty:
                    continue

                # normalize index name (get_var uses Capitalize of AMPL set names)
                phase_col = None
                for c in df.index.names:
                    if 'phase' in c.lower():
                        phase_col = c
                        break
                if phase_col is None:
                    self.results[var_name] = df
                    continue

                # filter to phases_wnd
                df = df.loc[
                    df.index.get_level_values(phase_col).isin(phases_wnd), :
                ]

                prev = self.results.get(var_name)
                if prev is None or (isinstance(prev, pd.DataFrame) and prev.empty):
                    self.results[var_name] = df
                else:
                    combined = pd.concat([prev, df])
                    self.results[var_name] = combined.loc[
                        ~combined.index.duplicated(keep='last')
                    ].sort_index()

                logging.info(f'  [collector] {var_name} accumulated: '
                             f'{len(self.results[var_name])} rows')

            except Exception as e:
                logging.warning(f'  [collector] Could not extract {var_name}: {e}')

        # ------------------------------------------------------------------
        # 3. Compute 2015_2021 phase costs (first window only).
        #    AMPL skips this phase because it is not in PHASE_WND.
        # ------------------------------------------------------------------
        if is_first_window and '2015_2021' in phases_wnd:
            self._compute_initial_phase_costs_from_ampl(pathway_model)

    # =========================================================================
    # COMPUTE INITIAL PHASE (2015_2021) COSTS FROM AMPL DATA
    # =========================================================================

    def _compute_initial_phase_costs_from_ampl(self, pathway_model) -> None:
        """
        Manually compute C_inv_phase_tech, C_op_phase_tech, C_op_phase_res
        for the initial phase '2015_2021' from available AMPL data.

        In perfect foresight, AMPL computes these because '2015_2021' is part
        of PHASE_WND.  In rolling horizon window 0, '2015_2021' is NOT in
        PHASE_WND, so the constraints that compute the phase-cost variables
        are never triggered for this phase.

        The formulas (from the AMPL model) are:

            C_inv_phase_tech[p,c,i] = F_new[p,c,i]
                * annualised_factor[p]
                * (c_inv[y_start,c,i] + c_inv[y_stop,c,i]) / 2

            C_op_phase_tech[p,c,i] = t_phase[p]
                * (C_maint[y_start,c,i] + C_maint[y_stop,c,i]) / 2
                * annualised_factor[p]

            C_op_phase_res[p,c,i] = t_phase[p]
                * (C_op[y_start,c,i] + C_op[y_stop,c,i]) / 2
                * annualised_factor[p]

        For phase 2015_2021: y_start = YEAR_2015, y_stop = YEAR_2021.

        Parameters
        ----------
        pathway_model : EsmcPathway
            The pathway model with a live AMPL session.
        """
        import numpy as np

        ampl = pathway_model.esom.ampl
        PHASE_NAME = '2015_2021'
        Y_START = 'YEAR_2015'
        Y_STOP  = 'YEAR_2021'
        T_PHASE = 6  # 2021 - 2015

        logging.info(f'  [collector] Computing initial phase ({PHASE_NAME}) '
                     f'costs from AMPL data ...')

        # ---- Read i_rate and compute annualised_factor -------------------
        try:
            i_rate = float(ampl.getValue('i_rate'))
        except Exception:
            i_rate = 0.015  # AMPL default
            logging.warning(f'  [collector] Could not read i_rate, '
                            f'using default {i_rate}')

        # diff_2015_phase is the number of years from 2015 to the
        # mid-point of the phase.  For 2015_2021: (2015+2021)/2 - 2015 = 3.
        # We try to read it from AMPL first; if unavailable, compute.
        try:
            diff_2015 = float(
                ampl.getParameter('diff_2015_phase')
                    .getValues()
                    .toDict()
                    .get(PHASE_NAME, 3.0)
            )
        except Exception:
            diff_2015 = (2015 + 2021) / 2 - 2015  # = 3.0

        annualised_factor = 1.0 / (1.0 + i_rate) ** diff_2015
        logging.info(f'  [collector]   i_rate={i_rate}, '
                     f'diff_2015={diff_2015}, '
                     f'annualised_factor={annualised_factor:.6f}')

        # ==================================================================
        # C_inv_phase_tech
        # ==================================================================
        try:
            # F_new for 2015_2021 (declared for PHASE union {"2015_2021"})
            f_new_all = ampl.getVariable('F_new').getValues().toPandas()
            val_col = [c for c in f_new_all.columns if '.val' in str(c)][0]
            f_new_all = f_new_all[[val_col]].rename(columns={val_col: 'F_new'})

            n_lvl = f_new_all.index.nlevels
            is_multi = (n_lvl == 3)  # (Phase, Region, Tech)
            idx_names = (['Phase', 'Regions', 'Technologies'] if is_multi
                         else ['Phase', 'Technologies'])
            f_new_all.index.names = idx_names

            # Filter to 2015_2021
            f_new_init = f_new_all.loc[
                f_new_all.index.get_level_values('Phase') == PHASE_NAME
            ].copy()

            if f_new_init.empty:
                logging.warning('  [collector]   F_new for 2015_2021 is empty')
            else:
                # c_inv parameter  {YEARS, REGIONS, TECHNOLOGIES}
                c_inv_param = pathway_model.esom.get_param('c_inv')

                # Extract c_inv for y_start and y_stop
                if 'Years' in c_inv_param.index.names:
                    yr_level = 'Years'
                elif 'YEARS' in c_inv_param.index.names:
                    yr_level = 'YEARS'
                else:
                    yr_level = c_inv_param.index.names[0]

                c_inv_start = c_inv_param.loc[
                    c_inv_param.index.get_level_values(yr_level) == Y_START
                ].copy()
                c_inv_stop = c_inv_param.loc[
                    c_inv_param.index.get_level_values(yr_level) == Y_STOP
                ].copy()

                # Drop year level to align with F_new index
                c_inv_start = c_inv_start.droplevel(yr_level)
                c_inv_stop  = c_inv_stop.droplevel(yr_level)

                # Get the value column name
                c_inv_col = c_inv_start.columns[0]

                # Align: drop Phase from f_new_init to match c_inv index
                f_new_vals = f_new_init.droplevel('Phase')

                # Compute: F_new * annualised_factor * (c_inv_start + c_inv_stop) / 2
                c_inv_avg = (c_inv_start[c_inv_col].reindex(f_new_vals.index, fill_value=0)
                             + c_inv_stop[c_inv_col].reindex(f_new_vals.index, fill_value=0)) / 2.0

                c_inv_phase = f_new_vals['F_new'] * annualised_factor * c_inv_avg
                c_inv_phase = c_inv_phase.dropna()
                c_inv_phase = c_inv_phase[c_inv_phase.abs() > 1e-6]

                if len(c_inv_phase) > 0:
                    # Rebuild MultiIndex with Phase
                    new_idx = pd.MultiIndex.from_arrays(
                        [[PHASE_NAME] * len(c_inv_phase)] +
                        [c_inv_phase.index.get_level_values(n)
                         for n in c_inv_phase.index.names],
                        names=['Phase'] + list(c_inv_phase.index.names)
                    )
                    c_inv_df = pd.DataFrame(
                        {'C_inv_phase_tech': c_inv_phase.values},
                        index=new_idx
                    )

                    # Accumulate
                    prev = self.results.get('C_inv_phase_tech')
                    if prev is None or (isinstance(prev, pd.DataFrame) and prev.empty):
                        self.results['C_inv_phase_tech'] = c_inv_df
                    else:
                        combined = pd.concat([c_inv_df, prev])
                        self.results['C_inv_phase_tech'] = combined.loc[
                            ~combined.index.duplicated(keep='last')
                        ].sort_index()

                    logging.info(f'  [collector]   C_inv_phase_tech for '
                                 f'{PHASE_NAME}: {len(c_inv_df)} rows, '
                                 f'total={c_inv_df["C_inv_phase_tech"].sum():.2f}')

        except Exception as e:
            logging.warning(f'  [collector] Could not compute C_inv_phase_tech '
                            f'for {PHASE_NAME}: {e}')

        # ==================================================================
        # C_op_phase_tech  (t_phase * avg(C_maint) * annualised_factor)
        # ==================================================================
        try:
            c_maint_var = pathway_model.esom.get_var('C_maint')
            if not c_maint_var.empty:
                yr_level = ([n for n in c_maint_var.index.names
                             if 'year' in n.lower()] or [c_maint_var.index.names[0]])[0]

                c_maint_start = c_maint_var.loc[
                    c_maint_var.index.get_level_values(yr_level) == Y_START
                ].copy()
                c_maint_stop = c_maint_var.loc[
                    c_maint_var.index.get_level_values(yr_level) == Y_STOP
                ].copy()

                c_maint_start = c_maint_start.droplevel(yr_level)
                c_maint_stop  = c_maint_stop.droplevel(yr_level)

                maint_col = c_maint_start.columns[0]

                # Compute average; start might be 0 for base year
                avg_maint = (c_maint_start[maint_col].reindex(
                                 c_maint_stop.index, fill_value=0)
                             + c_maint_stop[maint_col]) / 2.0

                c_op_phase = T_PHASE * avg_maint * annualised_factor
                c_op_phase = c_op_phase.dropna()
                c_op_phase = c_op_phase[c_op_phase.abs() > 1e-6]

                if len(c_op_phase) > 0:
                    new_idx = pd.MultiIndex.from_arrays(
                        [[PHASE_NAME] * len(c_op_phase)] +
                        [c_op_phase.index.get_level_values(n)
                         for n in c_op_phase.index.names],
                        names=['Phase'] + list(c_op_phase.index.names)
                    )
                    c_op_tech_df = pd.DataFrame(
                        {'C_op_phase_tech': c_op_phase.values},
                        index=new_idx
                    )

                    prev = self.results.get('C_op_phase_tech')
                    if prev is None or (isinstance(prev, pd.DataFrame) and prev.empty):
                        self.results['C_op_phase_tech'] = c_op_tech_df
                    else:
                        combined = pd.concat([c_op_tech_df, prev])
                        self.results['C_op_phase_tech'] = combined.loc[
                            ~combined.index.duplicated(keep='last')
                        ].sort_index()

                    logging.info(f'  [collector]   C_op_phase_tech for '
                                 f'{PHASE_NAME}: {len(c_op_tech_df)} rows')

        except Exception as e:
            logging.warning(f'  [collector] Could not compute C_op_phase_tech '
                            f'for {PHASE_NAME}: {e}')

        # ==================================================================
        # C_op_phase_res  (t_phase * avg(C_op) * annualised_factor)
        # ==================================================================
        try:
            c_op_var = pathway_model.esom.get_var('C_op')
            if not c_op_var.empty:
                yr_level = ([n for n in c_op_var.index.names
                             if 'year' in n.lower()] or [c_op_var.index.names[0]])[0]

                c_op_start = c_op_var.loc[
                    c_op_var.index.get_level_values(yr_level) == Y_START
                ].copy()
                c_op_stop = c_op_var.loc[
                    c_op_var.index.get_level_values(yr_level) == Y_STOP
                ].copy()

                c_op_start = c_op_start.droplevel(yr_level)
                c_op_stop  = c_op_stop.droplevel(yr_level)

                op_col = c_op_start.columns[0]

                avg_op = (c_op_start[op_col].reindex(
                              c_op_stop.index, fill_value=0)
                          + c_op_stop[op_col]) / 2.0

                c_op_res_phase = T_PHASE * avg_op * annualised_factor
                c_op_res_phase = c_op_res_phase.dropna()
                c_op_res_phase = c_op_res_phase[c_op_res_phase.abs() > 1e-6]

                if len(c_op_res_phase) > 0:
                    new_idx = pd.MultiIndex.from_arrays(
                        [[PHASE_NAME] * len(c_op_res_phase)] +
                        [c_op_res_phase.index.get_level_values(n)
                         for n in c_op_res_phase.index.names],
                        names=['Phase'] + list(c_op_res_phase.index.names)
                    )
                    c_op_res_df = pd.DataFrame(
                        {'C_op_phase_res': c_op_res_phase.values},
                        index=new_idx
                    )

                    prev = self.results.get('C_op_phase_res')
                    if prev is None or (isinstance(prev, pd.DataFrame) and prev.empty):
                        self.results['C_op_phase_res'] = c_op_res_df
                    else:
                        combined = pd.concat([c_op_res_df, prev])
                        self.results['C_op_phase_res'] = combined.loc[
                            ~combined.index.duplicated(keep='last')
                        ].sort_index()

                    logging.info(f'  [collector]   C_op_phase_res for '
                                 f'{PHASE_NAME}: {len(c_op_res_df)} rows')

        except Exception as e:
            logging.warning(f'  [collector] Could not compute C_op_phase_res '
                            f'for {PHASE_NAME}: {e}')

    # =========================================================================
    # RESTORE TO PATHWAY MODEL  (call before post-processing / plots)
    # =========================================================================

    def restore_to_pathway_model(self, pathway_model) -> None:
        """
        Restore the full-pathway results into pathway_model so that all
        post-processing and plot methods work without a live AMPL session.

        What this does:
          1. Copies self.results -> pathway_model.results
          2. Installs a MockEsom on pathway_model.esom so that methods that
             call self.esom.get_var(name) receive the stored DataFrame.

        Call this after load_from_pkl() or after the full rolling-horizon loop,
        just before running collect_new_old_decom / graph_cost_inv_phase_tech /
        graph_cost_op_phase.

        Parameters
        ----------
        pathway_model : EsmcPathway
            The pathway model (does not need an active AMPL session).
        """

        # 1. Restore pathway_model.results
        pathway_model.results = {k: v for k, v in self.results.items()
                                 if v is not None}

        # 2. Install MockEsom so that self.esom.get_var() works
        stored = self.results   # captured by closure

        class _MockAmpl:
            """Minimal AMPL mock — only used by collect_new_old_decom (getVariable)."""

            def getValue(self, name):
                return 'cached'

            def getVariable(self, var_name):
                """
                Return a mock AMPL variable whose getValues().toPandas() yields
                the raw format expected by collect_new_old_decom.
                For F_new / F_old / F_decom we reconstruct from New_old_decom.
                """

                class _MockDF:
                    def __init__(self, df):
                        self._df = df
                    def toPandas(self):
                        return self._df

                class _MockVar:
                    def __init__(self, df):
                        self._df = df
                    def getValues(self):
                        return _MockDF(self._df)

                nod = stored.get('New_old_decom')
                if nod is None or nod.empty:
                    return _MockVar(pd.DataFrame())

                col_map = {'F_new': 'F_new.val',
                           'F_old': 'F_old.val',
                           'F_decom': 'F_decom.val'}
                if var_name in col_map and col_map[var_name].replace('.val', '') in nod.columns:
                    raw_col = col_map[var_name].replace('.val', '')
                    df_out = nod[[raw_col]].rename(columns={raw_col: col_map[var_name]})
                    return _MockVar(df_out)

                return _MockVar(pd.DataFrame())

            def getSet(self, set_name):
                class _MockSet:
                    def getValues(self_inner):
                        return []
                    def __iter__(self_inner):
                        return iter([])
                return _MockSet()

            def close(self):
                pass

        class _MockEsom:
            def __init__(self):
                self.ampl = _MockAmpl()

            def get_var(self, var_name):
                df = stored.get(var_name)
                if df is not None and not df.empty:
                    return df.copy()
                return pd.DataFrame()

            def get_param(self, param_name):
                return pd.DataFrame()

        pathway_model.esom = _MockEsom()
        logging.info('[collector] Full pathway results restored to pathway_model')
        logging.info(f'  Results keys: {[k for k, v in self.results.items() if v is not None]}')

    # =========================================================================
    # LOAD FROM PKL  (classmethod — for --skip-optimize)
    # =========================================================================

    @classmethod
    def load_from_pkl(cls, output_file) -> 'EsmcCollector':
        """
        Load a previously saved _Results.pkl and return a populated collector.

        This is used by the --skip-optimize path to restore all rolling-horizon
        results without re-running the optimization.

        Parameters
        ----------
        output_file : Path or str
            Path to the _Results.pkl file saved by save_results().

        Returns
        -------
        EsmcCollector
            A collector instance with results, results_all and metadata populated.

        Raises
        ------
        FileNotFoundError
            If the pkl file does not exist.
        """
        output_file = Path(output_file)
        if not output_file.exists():
            raise FileNotFoundError(
                f'Results file not found: {output_file}\n'
                f'Run the full optimization first (without --skip-optimize).')

        with open(output_file, 'rb') as f:
            data = pickle.load(f)

        # Create a minimal collector (no esmc_pre needed for post-processing)
        obj = cls.__new__(cls)
        obj.output_file   = output_file
        obj.pth_output_all = output_file.parent.parent
        obj.results       = data.get('results', {})
        obj.results_all   = data.get('results_all', {})
        obj.expl_text     = data.get('metadata', {}).get('expl_text', '')
        obj._windows_processed = data.get('metadata', {}).get('windows_processed', [])
        obj.esmc_pre      = None   # not needed for post-processing
        obj.hourly_results = {}

        meta = data.get('metadata', {})
        logging.info(f'[collector] Loaded results from {output_file}')
        logging.info(f'  Created  : {meta.get("timestamp", "unknown")}')
        logging.info(f'  Windows  : {meta.get("n_windows", "unknown")}')
        logging.info(f'  Results  : {[k for k, v in obj.results.items() if v is not None]}')

        return obj

    def print_summary(self) -> None:
        """Print a summary of collected results."""
        summary = self.get_summary()

        print(f"\n{'='*60}")
        print(f"Results Collection Summary")
        print(f"{'='*60}")
        print(f"Windows processed: {summary['windows_processed']}")
        print(f"Total rows collected: {summary['total_rows']}")
        print(f"\nResult types: {len(summary['result_types'])}")
        print(f"{'-'*60}")

        for k, info in summary['result_types'].items():
            print(f"  {k}:")
            print(f"    Rows: {info['rows']}, Columns: {info['columns']}")
            print(f"    Index: {info['index_names']}")

        print(f"{'='*60}")