#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simplified EsmcPathway class for multi-year optimization
Works with pre-generated .dat files instead of regenerating from scratch
"""

import logging
import pandas as pd
import numpy as np
from pathlib import Path
import re
import json
from esmc.utils.opti_probl import OptiProbl
from esmc.common import CSV_SEPARATOR, AMPL_SEPARATOR
# Additional imports for TD generation
from esmc.utils.region import Region
from esmc.preprocessing.temporal_aggregation import TemporalAggregation
import esmc.preprocessing.dat_print as dp
import copy
import shutil
import csv




class SimpleTemporalAggregation:
    """
    Simplified temporal aggregation class for pathway analysis
    Only includes the essential methods needed to expand TD data to yearly data
    """
    
    def __init__(self, dat_dir: Path, Nbr_TD=16):
        """
        Initialize with minimum required data
        
        Parameters
        ----------
        dat_dir : Path
            Directory containing TD_of_days file
        Nbr_TD : int
            Number of typical days
        """
        self.dat_dir = dat_dir
        self.Nbr_TD = Nbr_TD
        self.td_of_days = None
        self.t_h_td = None
        self.td_count = None
        
        # Try to read TD_of_days file
        try:
            self.td_of_days = self.read_td_of_days()
            self.generate_t_h_td()
            logging.info(f'SimpleTemporalAggregation initialized with {Nbr_TD} typical days')
        except Exception as e:
            logging.warning(f'Could not initialize temporal aggregation: {e}')
            logging.warning('Will use direct sum without temporal expansion')
    
    def read_td_of_days(self, td_file=None):
        """
        Reads the file containing the TD_of_days mapping
        
        Parameters
        ----------
        td_file : Path, optional
            Path to TD_of_days file. If None, uses default location.
            
        Returns
        -------
        pd.DataFrame
            DataFrame with TD_of_days for each day of the year
        """
        if td_file is None:
            # Try multiple possible locations
            possible_paths = [
                self.dat_dir / f'TD_of_days_{self.Nbr_TD}.out',
                self.dat_dir.parent / '00_td_dat' / f'TD_of_days_{self.Nbr_TD}.out',
            ]
            
            td_file = None
            for path in possible_paths:
                if path.exists():
                    td_file = path
                    break
            
            if td_file is None:
                raise FileNotFoundError(f'TD_of_days file not found in: {possible_paths}')
        
        logging.info(f'Reading typical days from {td_file}')
        td_of_days = pd.read_csv(td_file, header=None)
        td_of_days.columns = ['TD_of_days']
        td_of_days.index = np.arange(1, 366)
        return td_of_days
    
    def generate_t_h_td(self):
        """
        Generate t_h_td and td_count dataframes
        
        t_h_td contains the mapping from year hours to typical day hours
        td_count contains the number of days each TD represents
        """
        if self.td_of_days is None:
            raise ValueError('td_of_days must be loaded first')
        
        # Get td_of_days and add day numbers
        td_of_days = self.td_of_days.copy()
        td_of_days['day'] = np.arange(1, 366, 1)
        
        # Compute number of days represented by each TD
        td_count = td_of_days.groupby('TD_of_days').count()
        td_count = td_count.reset_index().rename(columns={'day': '#days'})
        td_count['TD_number'] = np.arange(1, self.Nbr_TD + 1)
        
        # Build T_H_TD matrix (maps year hours to TD hours)
        t_h_td = pd.DataFrame(
            np.repeat(td_of_days['TD_of_days'].values, 24, axis=0),
            columns=['TD_of_days']
        )
        
        # Create mapping from TD_of_days to TD_number
        map_td = dict(zip(td_count['TD_of_days'], np.arange(1, self.Nbr_TD + 1)))
        t_h_td['TD_number'] = t_h_td['TD_of_days'].map(map_td)
        t_h_td['H_of_D'] = np.resize(np.arange(1, 25), t_h_td.shape[0])
        t_h_td['H_of_Y'] = np.arange(1, 8761)
        
        self.t_h_td = t_h_td
        self.td_count = td_count
        
        logging.info('t_h_td and td_count generated')
        return
    
    def from_td_to_year(self, ts_td):
        """
        Converts time series on TDs to yearly time series
        
        This expands data from typical days (e.g., 16 TDs * 24 hours = 384 rows)
        to full year (365 days * 24 hours = 8760 rows)
        
        Parameters
        ----------
        ts_td : pd.DataFrame
            Multiindex dataframe of hourly data for each hour of each TD.
            The index should be (TD_number, H_of_D) or (Typical_days, Hours)
            
        Returns
        -------
        pd.DataFrame
            Yearly time series with 8760 rows
        """
        if self.t_h_td is None:
            self.generate_t_h_td()
        
        # Get the mapping columns
        td_h = self.t_h_td.loc[:, ['TD_number', 'H_of_D']]
        
        # The index might be named differently, so we need to check
        index_names = ts_td.index.names
        
        # Rename to match what we need
        ts_td_copy = ts_td.reset_index()
        
        # Find which columns correspond to TD and Hour
        td_col = None
        h_col = None
        
        for col in ts_td_copy.columns:
            col_lower = str(col).lower()
            if 'typical' in col_lower or 'td' in col_lower:
                td_col = col
            elif 'hour' in col_lower or 'h_of' in col_lower:
                h_col = col
        
        if td_col is None or h_col is None:
            # Try by position if names don't match
            if len(index_names) >= 2:
                td_col = ts_td_copy.columns[0]
                h_col = ts_td_copy.columns[1]
        
        # Rename to standard names for merge
        rename_dict = {td_col: 'TD_number', h_col: 'H_of_D'}
        ts_td_copy = ts_td_copy.rename(columns=rename_dict)
        
        # Merge to expand from TDs to full year
        ts_yr = td_h.merge(
            ts_td_copy,
            left_on=['TD_number', 'H_of_D'],
            right_on=['TD_number', 'H_of_D']
        ).sort_index()
        
        # Drop the TD and Hour columns, keep only the data columns
        ts_yr = ts_yr.drop(columns=['TD_number', 'H_of_D'])
        
        return ts_yr



class EsmcPathway:
    """
    Simplified pathway optimization class
    Reads existing .dat files and combines them with year indexing
    """
    
    def __init__(self, config, years_config, nbr_td=16):
        """
        Initialize pathway optimization framework
        
        This simplified version works with pre-generated .dat files
        rather than regenerating everything from scratch.
        """
        # Basic configuration
        self.regions_names = config['regions_names']
        self.case_study = config['case_study']
        self.gwp_limit_overall = config.get('gwp_limit_overall', True)
        self.re_share_primary = config.get('re_share_primary', None)
        self.f_perc = config.get('f_perc', True)
        
        # Pathway configuration
        self.years = years_config['years']
        self.phases = years_config['phases']
        self.pathway_params = years_config.get('pathway_params', {})
        self.nbr_td = nbr_td
        
        # Directory structure
        self.project_dir = Path('/home/pjimenez/ESMC_PTH/EnergyScope_Unif')
        self.data_dir = self.project_dir / 'Data'
        self.pathway_dir = self.data_dir / 'ALL_YEARS'
        
        # Output directory
        self.space_id = config.get('space_id', 'default')
        self.cs_dir = self.project_dir / 'case_studies' / self.space_id / f"PATHWAY_{self.case_study}"
        self.cs_dir.mkdir(parents=True, exist_ok=True)

        # Create empty dictionaries to be filled with main results
        self.results = dict.fromkeys(['TotalCost', 'Cost_breakdown', 'Gwp_breakdown', 
                                    'Transfer_capacity', 'Exchanges_year', 'Resources', 
                                    'Exch_freight_border', 'Exch_freight',
                                    'Assets', 'Sto_assets', 'Year_balance', 'Curt'])
        
        self.hourly_results = dict()
        
        self.results_all = dict.fromkeys(['TotalCost', 'Cost_breakdown', 'Gwp_breakdown', 
                                        'Transfer_capacity', 'Exchanges_year', 'Resources', 
                                        'Exch_freight', 'Assets', 'Sto_assets', 
                                        'Year_balance', 'Curt'])
        
        # Simplified temporal aggregation will be initialized after TD generation
        # This avoids errors when TD files don't exist yet
        self.ta = None
        
        logging.info(f'EsmcPathway initialized for years: {self.years}')
        logging.info(f'Output directory: {self.cs_dir}')
    def read_data_indep(self, year='2021'):
        """
        Read independent data from the specified year
        This data is used for TD generation
        
        Parameters
        ----------
        year : str
            Year to read data from (default: '2021')
        """
        logging.info(f'Reading independent data from year {year}')
        data_dir = self.data_dir / str(year)
        
        # Initialize data_indep dictionary with required keys
        self.data_indep = dict.fromkeys(['Misc_indep', 'Layers_in_out', 'Resources_indep',
                                         'Storage_characteristics', 'Storage_eff_in', 
                                         'Storage_eff_out'])
        
        self.data_indep['Misc_indep'] = self._get_default_time_series_mapping()
        
        logging.info('Independent data read successfully')
        return
    
    def _get_default_time_series_mapping(self):
        """Get default time series mapping
        
        IMPORTANT: 
        - res_params: Only for 1-to-1 where TIME SERIES NAME == TECHNOLOGY NAME
        - res_mult_params: For 1-to-many or when names differ
        """
        return {
            'time_series_mapping': {
                'eud_params': {
                    'ELECTRICITY': 'electricity_time_series',
                    'HEAT_LOW_T_SH': 'heating_time_series',
                    'SPACE_COOLING': 'cooling_time_series',
                    'MOBILITY_PASSENGER': 'mob_pass_time_series',
                    'MOBILITY_FREIGHT': 'mob_freight_time_series'
                },
                'res_params': {
                    # Only include if tech name EQUALS time series name
                    'WIND_ONSHORE': 'WIND_ONSHORE',
                    'WIND_OFFSHORE': 'WIND_OFFSHORE',
                    'HYDRO_DAM': 'HYDRO_DAM',
                    'HYDRO_RIVER': 'HYDRO_RIVER'
                },
                'res_mult_params': {
                    # Multiple techs use same time series OR names differ
                    'PV': ['PV_ROOFTOP', 'PV_UTILITY'],  # PV time series -> 2 techs
                    'TIDAL': ['TIDAL_STREAM', 'TIDAL_RANGE'],
                    'SOLAR': ['DHN_SOLAR', 'DEC_SOLAR', 'PT_COLLECTOR', 
                             'ST_COLLECTOR', 'STIRLING_DISH']
                }
            }
        }
    
    def _validate_time_series_mapping(self):
        """
        Validate time_series_mapping against available technologies
        
        CRITICAL: For res_params, the time series name (KEY) must ALSO exist 
        as a technology name because compute_cell_w uses it to index Technologies.
        
        Returns validated mapping with only existing technologies.
        """
        if not hasattr(self, 'regions') or not self.regions:
            logging.warning('No regions initialized, using full time_series_mapping')
            return self.data_indep['Misc_indep']['time_series_mapping']
        
        # Get the full mapping
        full_mapping = self.data_indep['Misc_indep']['time_series_mapping'].copy()
        
        # Get available technologies from reference region
        ref_region = list(self.regions.values())[0]
        if ref_region is None:
            logging.warning('Reference region not initialized, using full time_series_mapping')
            return full_mapping
        
        available_techs = list(ref_region.data['Technologies'].index)
        logging.info(f'Available technologies: {len(available_techs)} technologies found')
        
        # Validate res_params
        # KEY (time series name) must ALSO be a technology name (same name required)
        validated_res_params = {}
        for ts_name, tech_name in full_mapping['res_params'].items():
            # Check if the KEY (time series name) exists as a technology
            if ts_name in available_techs:
                validated_res_params[ts_name] = tech_name
                logging.info(f'  res_params: {ts_name} -> {tech_name} [OK]')
            else:
                logging.warning(f'  res_params: {ts_name} not found as technology, moving to res_mult_params')
                # Try to find technologies that might use this time series
                # Look for technologies that start with the time series name
                matching_techs = [t for t in available_techs if ts_name in t]
                if matching_techs:
                    logging.info(f'    Found matching technologies: {matching_techs}')
                    # Add to res_mult_params instead
                    if 'res_mult_params' not in full_mapping:
                        full_mapping['res_mult_params'] = {}
                    full_mapping['res_mult_params'][ts_name] = matching_techs
        
        # Validate res_mult_params
        # Each technology in the VALUE list must exist
        validated_res_mult_params = {}
        for ts_name, tech_list in full_mapping['res_mult_params'].items():
            # Filter the list to only include existing technologies
            existing_techs = [t for t in tech_list if t in available_techs]
            if existing_techs:
                validated_res_mult_params[ts_name] = existing_techs
                logging.info(f'  res_mult_params: {ts_name} -> {existing_techs} [{len(existing_techs)}/{len(tech_list)} techs found]')
            else:
                logging.warning(f'  res_mult_params: No technologies found for {ts_name}, skipping')
        
        # Keep eud_params as is (demands should always exist)
        validated_mapping = {
            'eud_params': full_mapping['eud_params'],
            'res_params': validated_res_params,
            'res_mult_params': validated_res_mult_params
        }
        
        logging.info(f'Validation complete: {len(validated_res_params)} res_params, {len(validated_res_mult_params)} res_mult_params')
        
        return validated_mapping

    def init_regions(self, year='2021'):
        """
        Initialize regions from the specified year
        This is needed for TD generation
        
        Parameters
        ----------
        year : str
            Year to initialize regions from (default: '2021')
        """
        logging.info(f'Initialising regions from year {year}: ' + ', '.join(self.regions_names))
        data_dir = self.data_dir / str(year)
        
        # Create reference region
        ref_region_name = '02_REF_REGION'
        self.ref_region = Region(nuts=ref_region_name, data_dir=data_dir, ref_region=True)
        
        # Create regions dictionary
        self.regions = dict.fromkeys(self.regions_names, None)
        for r in self.regions_names:
            if r != ref_region_name:
                self.regions[r] = copy.deepcopy(self.ref_region)
                self.regions[r].__init__(nuts=r, data_dir=data_dir, ref_region=False)
            else:
                self.regions[r] = self.ref_region
        
        logging.info('Regions initialized successfully')
        return
    
    def init_ta(self, algo='kmedoid', ampl_path=None):
        """
        Initialize the temporal aggregator to generate TD files
        
        Parameters
        ----------
        algo : str
            Algorithm to use: 'kmedoid' to generate new TDs, 'read' to read existing
        ampl_path : str or None
            Path to AMPL executable
        """
        
        logging.info(f'Initializing TemporalAggregation with {algo} algorithm')
        
        # Auto-detect if TD files exist when algo='read'
        dat_dir = self.project_dir / 'case_studies' / self.space_id / '00_td_dat'
        dat_dir.mkdir(parents=True, exist_ok=True)
        
        td_file = dat_dir / f'TD_of_days_{self.nbr_td}.out'
        
        if algo == 'read' and not td_file.exists():
            logging.warning(f'TD file not found: {td_file}')
            logging.warning(f'Automatically switching to algo=kmedoid to generate TD files')
            algo = 'kmedoid'
        
        if algo == 'kmedoid':
            logging.info(f'Will generate TD files using kmedoid clustering')
        else:
            logging.info(f'Will read existing TD files from {td_file}')
        
        # Validate and filter time_series_mapping based on available technologies
        validated_mapping = self._validate_time_series_mapping()
        
        # Initialize full temporal aggregation object with validated mapping
        self.ta_full = TemporalAggregation(
            self.regions, 
            dat_dir, 
            Nbr_TD=self.nbr_td, 
            algo=algo,
            ampl_path=ampl_path,
            time_series_mapping=validated_mapping
        )
        
        logging.info(f'Temporal aggregation initialized with {self.nbr_td} typical days')
        logging.info(f'Time series error: {self.ta_full.e_ts}')
        
        # Initialize SimpleTemporalAggregation for result extraction
        # (after TD files have been generated by TemporalAggregation above)
        try:
            self.ta = SimpleTemporalAggregation(dat_dir=dat_dir, Nbr_TD=self.nbr_td)
            logging.info('SimpleTemporalAggregation initialized successfully')
        except Exception as e:
            logging.warning(f'Could not initialize SimpleTemporalAggregation: {e}')
            logging.warning('Will use ta_full for temporal expansion')
            self.ta = None
        
        return
    
    def print_td_data(self):
        """
        Print typical day data to reg_{nbr_td}TD.dat file
        This generates the TD file needed for the pathway optimization
        """
        # File to print to (in PATHWAY folder)
        # Path: /home/pjimenez/ESMC_PTH/EnergyScope_Unif/case_studies/BO/PATHWAY_NZE_2015_2050/reg_16TD.dat
        dat_file = self.cs_dir / f'reg_{self.nbr_td}TD.dat'
        
        logging.info(f'Printing TD data to {dat_file}')
        
        # PRELIMINARY COMPUTATIONS
        # Generate t_h_td and td_count, and compute rescaled typical days ts and peak factors
        self.ta_full.generate_t_h_td()
        
        peak_sh_factor = pd.DataFrame(0.0, index=self.regions_names, columns=['peak_sh_factor'])
        peak_sc_factor = pd.DataFrame(0.0, index=self.regions_names, columns=['peak_sc_factor'])
        
        for r in self.regions:
            self.regions[r].rescale_td_ts(self.ta_full.td_count)
            self.regions[r].compute_peak_sh_and_sc()
            peak_sh_factor.loc[r, 'peak_sh_factor'] = self.regions[r].peak_sh_factor
            peak_sc_factor.loc[r, 'peak_sc_factor'] = self.regions[r].peak_sc_factor
        
        # Prepare t_h_td for printing
        t_h_td = self.ta_full.t_h_td.copy()
        t_h_td['par_l'] = '('
        t_h_td['par_r'] = ')'
        t_h_td['comma1'] = ','
        t_h_td['comma2'] = ','
        t_h_td = t_h_td[['par_l', 'H_of_Y', 'comma1', 'H_of_D', 'comma2', 'TD_number', 'par_r']]
        
        # PRINTING
        # Print header
        header_file = self.project_dir / 'esmc' / 'energy_model' / 'headers' / 'header_td_data.txt'
        dp.print_header(dat_file=dat_file, header_file=header_file)
        
        # Print set T_H_TD
        dp.newline(dat_file, ['set T_H_TD := 		'])
        t_h_td.to_csv(dat_file, sep=AMPL_SEPARATOR, header=False, index=False, 
                     mode='a', quoting=csv.QUOTE_NONE)
        dp.end_table(dat_file)
        
        # Print parameters depending on TD
        dp.newline(dat_file, ['# -----------------------------', 
                              '# PARAMETERS DEPENDING ON NUMBER OF TYPICAL DAYS : ',
                              '# -----------------------------', ''])
        
        # Print nbr_tds
        dp.print_param(param=self.nbr_td, out_path=dat_file, name='nbr_tds')
        
        # Print peak factors
        dp.print_df(df=dp.ampl_syntax(peak_sh_factor), out_path=dat_file, name='param ')
        dp.print_df(df=dp.ampl_syntax(peak_sc_factor), out_path=dat_file, name='param ')
        
        # Get time series mappings
        eud_params = self.ta_full.time_series_mapping['eud_params']
        res_params = self.ta_full.time_series_mapping['res_params']
        res_mult_params = self.ta_full.time_series_mapping['res_mult_params']
        
        # Print EUD timeseries param
        for i in eud_params.keys():
            dp.newline(out_path=dat_file, comment=['param ' + eud_params[i] + ' :='])
            for r in self.regions:
                dp.print_df(
                    df=dp.ampl_syntax(self.regions[r].ts_td.loc[(i, slice(None)), :].droplevel(level=0)),
                    out_path=dat_file,
                    name='["' + r + '",*,*] : ', 
                    end_table=False
                )
            dp.end_table(out_path=dat_file)
        
        # Print c_p_t param
        dp.newline(out_path=dat_file, comment=['param c_p_t:='])
        
        # Print c_p_t part where 1 ts => 1 tech
        for i in res_params.keys():
            for r in self.regions:
                dp.print_df(
                    df=dp.ampl_syntax(self.regions[r].ts_td.loc[(i, slice(None)), :].droplevel(level=0)),
                    out_path=dat_file,
                    name='["' + res_params[i] + '","' + r + '",*,*] :',
                    end_table=False
                )
        
        # Print c_p_t part where 1 ts => more than 1 tech
        for i in res_mult_params.keys():
            for j in res_mult_params[i]:
                for r in self.regions:
                    dp.print_df(
                        df=dp.ampl_syntax(self.regions[r].ts_td.loc[(i, slice(None)), :].droplevel(level=0)),
                        out_path=dat_file, 
                        name='["' + j + '","' + r + '",*,*] :',
                        end_table=False
                    )
        
        dp.end_table(out_path=dat_file)
        
        logging.info(f'TD data written to {dat_file}')
        
        # Mark that TD was generated (to skip copying later)
        self._td_generated = True
        
        # Initialize SimpleTemporalAggregation now that TD files exist
        try:
            dat_dir = self.project_dir / 'case_studies' / self.space_id / '00_td_dat'
            self.ta = SimpleTemporalAggregation(dat_dir=dat_dir, Nbr_TD=self.nbr_td)
            logging.info('SimpleTemporalAggregation initialized for result extraction')
        except Exception as e:
            logging.warning(f'Could not initialize SimpleTemporalAggregation: {e}')
            logging.warning('Will use ta_full.from_td_to_year() for temporal expansion')
        
        
        return
    
    def process_pathway_data(self):
        """
        Main method to process all pathway data
        
        This reads existing .dat files from each year directory and combines them
        with appropriate year indexing for the pathway model.
        """
        logging.info('Processing pathway data from existing .dat files')
        
        # Note: TD files are generated via init_ta() and print_td_data(), not copied from ALL_YEARS
        
        # Process year-specific .dat files using the correct functions
        self._process_indep_dat()
        self.process_reg_demands_dat_for_pathway()
        self._process_resources_dat()
        self._process_technologies_dat()
        self._process_exchanges_dat()
        self._process_misc_dat()
        self._process_storage_dat()
        
        # Copy pathway-specific files from ALL_YEARS (PES files, etc.)
        self._copy_pathway_files()
        
        logging.info('All pathway data processed successfully')

    def _process_indep_dat(self):
        """
        Procesa archivos indep.dat de cada aÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â±o y los combina con indexaciÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â³n temporal.
        Maneja cada parÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡metro segÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Âºn su formato especÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â­fico.
        """
        output_file = self.cs_dir / 'indep_pathway.dat'
        
        try:
            with open(output_file, 'w') as out_f:
                out_f.write('# Independent data for pathway optimization\n')
                out_f.write('# Generated from year-specific indep.dat files\n\n')
                
                # Primero, copiar los SETS del primer aÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â±o (son constantes)
                first_year_file = self.data_dir / str(self.years[0]) / 'indep.dat'
                if first_year_file.exists():
                    with open(first_year_file, 'r') as f:
                        content = f.read()
                    
                    # Extraer y escribir todos los SETS
                    sets_match = re.findall(r'set\s+(\w+(?:\s*\[[^\]]+\])?)\s*:?=([^;]+);', content, re.DOTALL)
                    if sets_match:
                        out_f.write('# SETS (constant across years)\n')
                        for set_name, set_values in sets_match:
                            out_f.write(f'set {set_name} := {set_values};\n')
                        out_f.write('\n')
                
                # =======================
                # 1. FORMATO SIMPLE
                # =======================
                # gwp_limit_overall - formato: param gwp_limit_overall := valor;
                out_f.write('# gwp_limit_overall by year\n')
                out_f.write('param gwp_limit_overall :=\n')
                for year in self.years:
                    indep_file = self.data_dir / str(year) / 'indep.dat'
                    if indep_file.exists():
                        with open(indep_file, 'r') as f:
                            content = f.read()
                        
                        match = re.search(r'param\s+gwp_limit_overall\s*:?=\s*(Infinity|\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s*;', content)
                        if match:
                            out_f.write(f'  YEAR_{year} {match.group(1)}\n')
                out_f.write(';\n\n')
                
                # =======================
                # 2. FORMATO COLUMNA ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã‚Â¦ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡NICA
                # =======================
                # loss_network y vehicule_capacity
                single_col_params = ['loss_network', 'vehicule_capacity']
                
                for param in single_col_params:
                    out_f.write(f'# {param} by year\n')
                    out_f.write(f'param {param} :=\n')
                    
                    for year in self.years:
                        year_file = self.data_dir / str(year) / 'indep.dat'
                        if year_file.exists():
                            with open(year_file, 'r') as f:
                                content = f.read()
                            
                            # Buscar el parÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡metro con formato "param : param_name :="
                            pattern = rf'param\s*:\s*{param}\s*:=(.*?);'
                            match = re.search(pattern, content, re.DOTALL)
                            if match:
                                data = match.group(1).strip()
                                lines = data.split('\n')
                                
                                for line in lines:
                                    if line.strip() and not line.startswith('#'):
                                        parts = line.split()
                                        if len(parts) >= 2:
                                            out_f.write(f'  YEAR_{year} {parts[0]} {parts[1]}\n')
                    
                    out_f.write(';\n\n')
                
                # =======================
                # 3. FORMATO MULTI-COLUMNA
                # =======================
                # gwp_op_exterior, c_op_exterior, co2_net
                multi_col_params = ['gwp_op_exterior', 'c_op_exterior', 'co2_net']
                
                # Procesar cada parÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡metro por separado
                for param_name in multi_col_params:
                    out_f.write(f'# {param_name} by year\n')
                    out_f.write(f'param {param_name} :=\n')
                    
                    for year in self.years:
                        year_file = self.data_dir / str(year) / 'indep.dat'
                        if year_file.exists():
                            with open(year_file, 'r') as f:
                                year_content = f.read()
                            
                            # Buscar parÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡metros con mÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Âºltiples columnas
                            pattern = r'param\s*:\s*((?:gwp_op_exterior|c_op_exterior|co2_net)(?:\s+(?:gwp_op_exterior|c_op_exterior|co2_net))*)\s*:=(.*?);'
                            match = re.search(pattern, year_content, re.DOTALL)
                            
                            if match:
                                param_names = match.group(1).split()
                                data_lines = match.group(2).strip().split('\n')
                                
                                # Encontrar el ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â­ndice de la columna para este parÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡metro
                                if param_name in param_names:
                                    col_index = param_names.index(param_name)
                                    
                                    # Procesar cada lÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â­nea de datos
                                    for line in data_lines:
                                        parts = line.split()
                                        if len(parts) > col_index + 1:  # +1 porque la primera parte es el recurso
                                            resource = parts[0]
                                            value = parts[col_index + 1]
                                            out_f.write(f'  YEAR_{year} {resource} {value}\n')
                    
                    out_f.write(';\n\n')
                
                # =======================
                # 4. FORMATO DOBLE COLUMNA
                # =======================
                # storage_availability y storage_losses
                storage_params = ['storage_availability', 'storage_losses']
                
                for param in storage_params:
                    out_f.write(f'# {param} by year\n')
                    out_f.write(f'param {param} :=\n')
                    
                    for year in self.years:
                        year_file = self.data_dir / str(year) / 'indep.dat'
                        if year_file.exists():
                            with open(year_file, 'r') as f:
                                year_content = f.read()
                            
                            # Buscar el parÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡metro con formato "param : param1 param2 :="
                            pattern = rf'param\s*:\s*(?:\w+\s+)*{param}(?:\s+\w+)*\s*:=(.*?);'
                            match = re.search(pattern, year_content, re.DOTALL)
                            
                            if match:
                                # Buscar el encabezado para saber en quÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â© columna estÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡ el parÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡metro
                                header_match = re.search(
                                    rf'param\s*:\s*(.*?{param}.*?)\s*:=', 
                                    year_content
                                )
                                if header_match:
                                    headers = header_match.group(1).split()
                                    if param in headers:
                                        col_index = headers.index(param)
                                        
                                        data = match.group(1).strip()
                                        lines = data.split('\n')
                                        
                                        for line in lines:
                                            parts = line.split()
                                            if len(parts) > col_index:
                                                tech = parts[0]
                                                value = parts[col_index + 1]  # +1 porque parts[0] es la tecnologÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â­a
                                                out_f.write(f'  YEAR_{year} {tech} {value}\n')
                    
                    out_f.write(';\n\n')
                
                # =======================
                # 5. FORMATO MATRIZ
                # =======================
                # layers_in_out, storage_eff_in, storage_eff_out
                matrix_params = ['layers_in_out', 'storage_eff_in', 'storage_eff_out']
                
                for param in matrix_params:
                    out_f.write(f'# {param} matrix by year\n')
                    out_f.write(f'param {param} :=\n')
                    
                    for year in self.years:
                        year_file = self.data_dir / str(year) / 'indep.dat'
                        if year_file.exists():
                            with open(year_file, 'r') as f:
                                content = f.read()
                            
                            # Buscar el parÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡metro matriz
                            pattern = rf'param {param}\s*:\s*(.*?):=(.*?);'
                            match = re.search(pattern, content, re.DOTALL)
                            
                            if match:
                                headers = match.group(1).strip()
                                data = match.group(2).strip()
                                
                                if data:
                                    # Escribir con indexaciÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â³n por aÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â±o
                                    out_f.write(f'[YEAR_{year},*,*]:\n')
                                    out_f.write(f'  {headers} :=\n')
                                    
                                    for line in data.split('\n'):
                                        if line.strip() and not line.startswith('#'):
                                            out_f.write(f'  {line}\n')
                    
                    out_f.write(';\n\n')
                
            logging.info(f'Successfully processed indep.dat to {output_file}')
            
        except Exception as e:
            logging.error(f'Error processing indep.dat: {e}')
            import traceback
            traceback.print_exc()
    
    def process_reg_demands_dat_for_pathway(self):
        """
        Process regional demands data with year indexing.
        El parÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡metro end_uses_demand_year tiene formato matriz con 4 dimensiones:
        YEAR, REGION, END_USE_TYPE, SECTOR
        """
        output_file = self.cs_dir / 'reg_demands_pathway.dat'
        
        try:
            with open(output_file, 'w') as out_f:
                out_f.write('# Regional demands for pathway optimization\n')
                out_f.write('# Generated from year-specific reg_demands.dat files\n\n')
                
                # Primero escribir el set REGIONS (es constante para todos los aÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â±os)
                out_f.write(f'set REGIONS := {" ".join(self.regions_names)} ;\n\n')
                
                # Ahora procesar end_uses_demand_year
                out_f.write('param end_uses_demand_year :=\n')
                
                for year in self.years:
                    year_file = self.data_dir / str(year) / 'reg_demands.dat'
                    if not year_file.exists():
                        logging.warning(f'reg_demands.dat not found for year {year}, skipping')
                        continue
                    
                    logging.info(f'Processing reg_demands.dat for year {year}')
                    
                    with open(year_file, 'r') as f:
                        year_content = f.read()
                    
                    # Extraer el bloque completo de end_uses_demand_year
                    pattern = r'param end_uses_demand_year\s*:\s*(.*?):=(.*?)(?=;|\Z)'
                    match = re.search(pattern, year_content, re.DOTALL)
                    
                    if match:
                        # Extraer los nombres de los sectores (columnas)
                        sectors_line = match.group(1).strip()
                        sectors = sectors_line.split()
                        
                        # Extraer los datos
                        data_block = match.group(2).strip()
                        
                        # Procesar cada lÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â­nea de datos
                        for line in data_block.split('\n'):
                            line = line.strip()
                            if not line or line.startswith('#'):
                                continue
                            
                            parts = line.split()
                            if len(parts) >= 3:  # Al menos regiÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â³n, tipo de demanda y un valor
                                region = parts[0]
                                end_use_type = parts[1]
                                
                                # Solo procesar si es una regiÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â³n vÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡lida
                                if region in self.regions_names:
                                    # Los valores comienzan desde el ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â­ndice 2
                                    values = parts[2:]
                                    
                                    # Para cada sector y su valor correspondiente
                                    for i, sector in enumerate(sectors):
                                        if i < len(values):
                                            value = values[i]
                                            # Solo escribir si el valor es significativo
                                            try:
                                                float_value = float(value)
                                                if abs(float_value) > 1e-6:  # Evitar valores muy pequeÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â±os
                                                    out_f.write(f'  YEAR_{year} {region} {end_use_type} {sector} {value}\n')
                                            except ValueError:
                                                logging.warning(f'Could not parse value: {value} for {region} {end_use_type} {sector}')
                
                out_f.write(';\n\n')
                
                logging.info(f'Successfully processed reg_demands.dat to {output_file}')
                return True
                
        except Exception as e:
            logging.error(f'Error processing reg_demands.dat: {e}')
            import traceback
            traceback.print_exc()
            return False
        
    def _process_resources_dat(self):
        """
        Procesa archivos reg_resources.dat de cada aÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â±o y los combina con indexaciÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â³n temporal.
        Los parÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡metros avail_local, avail_exterior, gwp_op_local y c_op_local estÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡n en formato multi-columna.
        """
        output_file = self.cs_dir / 'reg_resources_pathway.dat'
        
        try:
            with open(output_file, 'w') as out_f:
                out_f.write('# Regional resources data for pathway optimization\n')
                out_f.write('# Generated from year-specific reg_resources.dat files\n\n')
                
                # Los cuatro parÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡metros que necesitamos procesar
                resource_params = ['avail_local', 'avail_exterior', 'gwp_op_local', 'c_op_local']
                
                # Procesar cada parÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡metro individualmente con indexaciÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â³n por aÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â±o
                for param_name in resource_params:
                    out_f.write(f'# {param_name} by year and region\n')
                    # CAMBIO CLAVE: Eliminamos "default 0" de la declaraciÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â³n
                    # Antes: out_f.write(f'param {param_name} default 0 :=\n')
                    # Ahora: simplemente declaramos el parÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡metro sin valor por defecto
                    out_f.write(f'param {param_name} :=\n')
                    
                    for year in self.years:
                        year_file = self.data_dir / str(year) / 'reg_resources.dat'
                        if year_file.exists():
                            with open(year_file, 'r') as f:
                                content = f.read()
                            
                            # Buscar la lÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â­nea con los nombres de parÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡metros
                            # Formato: param : avail_local avail_exterior gwp_op_local c_op_local :=
                            pattern = r'param\s*:\s*(.*?):='
                            header_match = re.search(pattern, content)
                            
                            if header_match:
                                headers = header_match.group(1).split()
                                
                                # Verificar que el parÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡metro estÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â© en los headers
                                if param_name in headers:
                                    param_index = headers.index(param_name)
                                    
                                    # Extraer los datos despuÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â©s del :=
                                    data_pattern = r':=(.*?);'
                                    data_match = re.search(data_pattern, content, re.DOTALL)
                                    
                                    if data_match:
                                        data_lines = data_match.group(1).strip().split('\n')
                                        
                                        for line in data_lines:
                                            if line.strip() and not line.startswith('#'):
                                                parts = line.split()
                                                # parts[0] = regiÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â³n, parts[1] = recurso, parts[2+] = valores
                                                if len(parts) >= param_index + 3:  # +2 por regiÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â³n y recurso, +1 para el ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â­ndice
                                                    region = parts[0]
                                                    resource = parts[1]
                                                    value = parts[param_index + 2]  # +2 porque ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â­ndices 0 y 1 son regiÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â³n y recurso
                                                    
                                                    # Escribir con formato: YEAR regiÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â³n recurso valor
                                                    out_f.write(f'  YEAR_{year} {region} {resource} {value}\n')
                    
                    out_f.write(';\n\n')
                
            logging.info(f'Successfully processed reg_resources.dat to {output_file}')
            
        except Exception as e:
            logging.error(f'Error processing reg_resources.dat: {e}')
            import traceback
            traceback.print_exc()
        
    def _process_technologies_dat(self):
        """
        Procesa archivos reg_technologies.dat de cada aÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â±o y los combina con indexaciÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â³n temporal.
        Todos los parÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡metros se escriben en un ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Âºnico archivo .dat
        """
        output_file = self.cs_dir / 'reg_technologies_pathway.dat'
        
        try:
            with open(output_file, 'w') as out_f:
                out_f.write('# Regional technologies data for pathway optimization\n')
                out_f.write('# Generated from year-specific reg_technologies.dat files\n\n')
                
                # Lista de parÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡metros en el orden que aparecen en el archivo
                tech_params = ['c_inv', 'c_maint', 'gwp_constr', 'lifetime', 'c_p', 
                            'fmin_perc', 'fmax_perc', 'f_min', 'f_max']
                
                # Para cada parÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡metro, procesar todos los aÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â±os
                for param_idx, param_name in enumerate(tech_params):
                    out_f.write(f'# {param_name} by year and region\n')
                    
                    # ELIMINADO: Ya no especificamos valores por defecto
                    # Simplemente declaramos el parÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡metro sin 'default'
                    out_f.write(f'param {param_name} :=\n')
                    
                    # Procesar cada aÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â±o
                    for year in self.years:
                        year_file = self.data_dir / str(year) / 'reg_technologies.dat'
                        if year_file.exists():
                            with open(year_file, 'r') as f:
                                content = f.read()
                            
                            # Buscar la tabla de parÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡metros
                            # El formato es: param : c_inv c_maint ... := datos
                            pattern = r'param\s*:\s*(c_inv\s+c_maint.*?):=(.*?);'
                            match = re.search(pattern, content, re.DOTALL)
                            
                            if match:
                                # Obtener los nombres de columnas
                                headers = match.group(1).split()
                                data_section = match.group(2).strip()
                                
                                # Verificar que el parÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡metro estÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡ en los headers
                                if param_name in headers:
                                    param_col_idx = headers.index(param_name)
                                    
                                    # Procesar cada lÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â­nea de datos
                                    for line in data_section.split('\n'):
                                        if line.strip() and not line.startswith('#'):
                                            parts = line.split()
                                            if len(parts) >= len(headers) + 2:  # regiÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â³n + tech + params
                                                region = parts[0]
                                                technology = parts[1]
                                                value = parts[2 + param_col_idx]  # +2 por regiÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â³n y tech
                                                
                                                # Manejar valores especiales como 'Infinity'
                                                if value == 'Infinity':
                                                    value = '1e6'  # Usar un valor grande en lugar de Infinity
                                                
                                                out_f.write(f'  YEAR_{year} {region} {technology} {value}\n')
                    
                    out_f.write(';\n\n')
                
            logging.info(f'Successfully processed technologies data to {output_file}')
            
        except Exception as e:
            logging.error(f'Error processing technologies data: {e}')
            import traceback
            traceback.print_exc()
        
    def _process_exchanges_dat(self):
        """
        Procesa archivos reg_exch.dat de cada aÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â±o y los combina con indexaciÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â³n temporal.
        Maneja cada parÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡metro segÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Âºn su formato especÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â­fico.
        """
        output_file = self.cs_dir / 'reg_exchanges_pathway.dat'
        
        try:
            with open(output_file, 'w') as out_f:
                out_f.write('# Regional exchanges for pathway optimization\n')
                out_f.write('# Generated from year-specific reg_exch.dat files\n\n')
                
                # =======================
                # 1. FORMATO SIMPLE
                # =======================
                # retro_gas_to_h2
                out_f.write('# retro_gas_to_h2 by year\n')
                out_f.write('param retro_gas_to_h2 :=\n')
                
                for year in self.years:
                    year_file = self.data_dir / str(year) / 'reg_exch.dat'
                    if year_file.exists():
                        with open(year_file, 'r') as f:
                            content = f.read()
                        
                        # Buscar retro_gas_to_h2
                        match = re.search(r'param\s+retro_gas_to_h2\s*:?=\s*(Infinity|[\d.]+(?:[eE][+-]?\d+)?)', content)
                        if match:
                            out_f.write(f'  YEAR_{year} {match.group(1)}\n')
                
                out_f.write(';\n\n')
                
                # =======================
                # 2. FORMATO DOBLE COLUMNA
                # =======================
                # dist - las distancias no cambian por aÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â±o, pero agregamos indexaciÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â³n para consistencia
                out_f.write('# Distance matrix by year\n')
                out_f.write('param dist :=\n')
                
                for year in self.years:
                    year_file = self.data_dir / str(year) / 'reg_exch.dat'
                    if year_file.exists():
                        with open(year_file, 'r') as f:
                            content = f.read()
                        
                        # Buscar dist
                        pattern = r'param\s*:\s*dist\s*:=(.*?);'
                        match = re.search(pattern, content, re.DOTALL)
                        if match:
                            data = match.group(1).strip()
                            lines = data.split('\n')
                            
                            for line in lines:
                                if line.strip() and not line.startswith('#'):
                                    parts = line.split()
                                    if len(parts) >= 3:  # regiÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â³n1 regiÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â³n2 distancia
                                        out_f.write(f'  YEAR_{year} {parts[0]} {parts[1]} {parts[2]}\n')
                
                out_f.write(';\n\n')
                
                # =======================
                # 3. FORMATO COLUMNA ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã‚Â¦ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡NICA
                # =======================
                # exchange_losses
                out_f.write('# Exchange losses by year\n')
                out_f.write('param exchange_losses :=\n')
                
                for year in self.years:
                    year_file = self.data_dir / str(year) / 'reg_exch.dat'
                    if year_file.exists():
                        with open(year_file, 'r') as f:
                            content = f.read()
                        
                        # Buscar exchange_losses
                        pattern = r'param\s*:\s*exchange_losses\s*:=(.*?);'
                        match = re.search(pattern, content, re.DOTALL)
                        if match:
                            data = match.group(1).strip()
                            lines = data.split('\n')
                            
                            for line in lines:
                                if line.strip() and not line.startswith('#'):
                                    parts = line.split()
                                    if len(parts) >= 2:  # recurso valor
                                        out_f.write(f'  YEAR_{year} {parts[0]} {parts[1]}\n')
                
                out_f.write(';\n\n')
                
                # lhv
                out_f.write('# Lower heating values by year\n')
                out_f.write('param lhv :=\n')
                
                for year in self.years:
                    year_file = self.data_dir / str(year) / 'reg_exch.dat'
                    if year_file.exists():
                        with open(year_file, 'r') as f:
                            content = f.read()
                        
                        # Buscar lhv
                        pattern = r'param\s*:\s*lhv\s*:=(.*?);'
                        match = re.search(pattern, content, re.DOTALL)
                        if match:
                            data = match.group(1).strip()
                            lines = data.split('\n')
                            
                            for line in lines:
                                if line.strip() and not line.startswith('#'):
                                    parts = line.split()
                                    if len(parts) >= 2:  # recurso valor
                                        out_f.write(f'  YEAR_{year} {parts[0]} {parts[1]}\n')
                
                out_f.write(';\n\n')
                
                # =======================
                # 4. FORMATO MULTI-COLUMNA COMPLEJO
                # =======================
                # tc_min y tc_max estÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡n en el mismo bloque
                tc_params = ['tc_min', 'tc_max']
                
                for param in tc_params:
                    out_f.write(f'# {param} by year\n')
                    # ELIMINADO: default 0
                    out_f.write(f'param {param} :=\n')
                    
                    for year in self.years:
                        year_file = self.data_dir / str(year) / 'reg_exch.dat'
                        if year_file.exists():
                            with open(year_file, 'r') as f:
                                content = f.read()
                            
                            # Buscar el bloque tc_min tc_max
                            pattern = r'param\s*:\s*tc_min\s+tc_max\s*:=(.*?);'
                            match = re.search(pattern, content, re.DOTALL)
                            
                            if match:
                                data = match.group(1).strip()
                                lines = data.split('\n')
                                
                                # Determinar quÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â© columna corresponde a este parÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡metro
                                col_index = 4 if param == 'tc_min' else 5  # tc_min estÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡ en col 4, tc_max en col 5
                                
                                for line in lines:
                                    if line.strip() and not line.startswith('#'):
                                        parts = line.split()
                                        # Formato: regiÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â³n1 regiÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â³n2 recurso network_type tc_min tc_max
                                        if len(parts) >= 6:
                                            value = parts[col_index]
                                            # Escribir con indexaciÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â³n por aÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â±o
                                            out_f.write(f'  YEAR_{year} {parts[0]} {parts[1]} {parts[2]} {parts[3]} {value}\n')
                    
                    out_f.write(';\n\n')
                
                # =======================
                # PROCESAR OTROS PARÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚ÂMETROS DEL ARCHIVO reg_exch.dat
                # =======================
                # Estos parÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡metros no estÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡n en el archivo de texto de ejemplo pero estÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡n en el cÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â³digo original
                # Los incluyo para mantener la funcionalidad completa
                
                # exch_capacity, exch_lifetime, exch_c_inv, exch_c_maint
                year_varying_params = ['exch_capacity', 'exch_lifetime', 'exch_c_inv', 'exch_c_maint']
                
                for param in year_varying_params:
                    # Verificar si el parÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡metro existe en algÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Âºn archivo
                    param_exists = False
                    for year in self.years:
                        year_file = self.data_dir / str(year) / 'reg_exch.dat'
                        if year_file.exists():
                            with open(year_file, 'r') as f:
                                if f'param {param}' in f.read():
                                    param_exists = True
                                    break
                    
                    if param_exists:
                        out_f.write(f'# {param} by year\n')
                        # ELIMINADO: default 0
                        out_f.write(f'param {param} :=\n')
                        
                        for year in self.years:
                            year_file = self.data_dir / str(year) / 'reg_exch.dat'
                            if year_file.exists():
                                with open(year_file, 'r') as f:
                                    content = f.read()
                                
                                # Buscar el parÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡metro
                                pattern = rf'param {param}.*?:=(.*?);'
                                match = re.search(pattern, content, re.DOTALL)
                                if match:
                                    data = match.group(1).strip()
                                    lines = data.split('\n')
                                    
                                    for line in lines:
                                        if line.strip() and not line.startswith('#'):
                                            # Agregar aÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â±o al principio
                                            out_f.write(f'  YEAR_{year} {line.strip()}\n')
                        
                        out_f.write(';\n\n')
                
            logging.info(f'Successfully processed reg_exch.dat to {output_file}')
            
        except Exception as e:
            logging.error(f'Error processing reg_exch.dat: {e}')
            import traceback
            traceback.print_exc()
        
    def _process_misc_dat(self):
        """
        Process miscellaneous regional data with proper year indexing
        """
        output_file = self.cs_dir / 'reg_misc_pathway.dat'
        
        try:
            with open(output_file, 'w') as out_f:
                out_f.write('# Miscellaneous regional data for pathway optimization\n\n')
                
                # 1. RWITHOUTDAM set (constante, no varÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â­a por aÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â±o)
                first_year_file = self.data_dir / str(self.years[0]) / 'reg_misc.dat'
                if first_year_file.exists():
                    with open(first_year_file, 'r') as f:
                        content = f.read()
                    
                    # Extraer RWITHOUTDAM
                    rwithoutdam_match = re.search(r'set\s+RWITHOUTDAM\s*:?=\s*([^;]+);', content, re.DOTALL)
                    if rwithoutdam_match:
                        out_f.write('set RWITHOUTDAM := ' + rwithoutdam_match.group(1) + ';\n\n')
                
                # 2. share_ned - FORMATO MATRIZ (regiÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â³n x productos)
                out_f.write('# share_ned by year\n')
                out_f.write('param share_ned :=\n')
                
                for year in self.years:
                    misc_file = self.data_dir / str(year) / 'reg_misc.dat'
                    if misc_file.exists():
                        with open(misc_file, 'r') as f:
                            content = f.read()
                        
                        # Buscar share_ned con su formato especÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â­fico
                        pattern = r'param share_ned\s*:\s*(.*?):=(.*?);'
                        match = re.search(pattern, content, re.DOTALL)
                        
                        if match:
                            headers = match.group(1).strip()  # HVC METHANOL AMMONIA
                            data = match.group(2).strip()
                            
                            if data:
                                # Escribir con indexaciÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â³n [YEAR,*,*]
                                out_f.write(f'[YEAR_{year},*,*]:\n')
                                out_f.write(f'  {headers} :=\n')
                                
                                for line in data.split('\n'):
                                    if line.strip() and not line.startswith('#'):
                                        out_f.write(f'  {line}\n')
                
                out_f.write(';\n\n')
                
                # 3. Otros parÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡metros - FORMATO MULTI-COLUMNA
                # Lista de todos los parÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡metros agrupados por tabla
                param_groups = [
                    ['gwp_limit', 'import_capacity', 're_share_primary', 'share_freight_boat_max'],
                    ['share_freight_boat_min', 'share_freight_road_max', 'share_freight_road_min', 'share_freight_train_max'],
                    ['share_freight_train_min', 'share_heat_dhn_max', 'share_heat_dhn_min', 'share_mobility_public_max'],
                    ['share_mobility_public_min', 'share_short_haul_flights_max', 'share_short_haul_flights_min', 'solar_area_ground'],
                    ['solar_area_ground_high_irr', 'solar_area_rooftop']
                ]
                
                # Procesar cada parÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡metro individualmente
                all_params = [param for group in param_groups for param in group]
                
                for param_name in all_params:
                    out_f.write(f'# {param_name} by year\n')
                    
                    # CAMBIO CRÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚ÂTICO: Eliminamos la especificaciÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â³n de valores por defecto
                    # Ya no distinguimos entre gwp_limit (que tenÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â­a default Infinity) y otros parÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡metros
                    # Los valores por defecto deben estar definidos en el archivo .mod
                    out_f.write(f'param {param_name} :=\n')
                    
                    for year in self.years:
                        misc_file = self.data_dir / str(year) / 'reg_misc.dat'
                        if misc_file.exists():
                            with open(misc_file, 'r') as f:
                                content = f.read()
                            
                            # Buscar en quÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â© tabla estÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡ este parÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡metro
                            param_found = False
                            for group in param_groups:
                                if param_name in group:
                                    # Crear el patrÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â³n para buscar esta tabla especÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â­fica
                                    params_pattern = '\\s+'.join(group)
                                    table_pattern = rf'param\s*:\s*{params_pattern}\s*:=(.*?);'
                                    
                                    match = re.search(table_pattern, content, re.DOTALL)
                                    if match:
                                        data = match.group(1).strip()
                                        lines = data.split('\n')
                                        
                                        # Encontrar el ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â­ndice de la columna para este parÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡metro
                                        col_index = group.index(param_name)
                                        
                                        # Procesar cada lÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â­nea
                                        for line in lines:
                                            if line.strip() and not line.startswith('#'):
                                                parts = line.split()
                                                if len(parts) > col_index + 1:  # +1 porque la primera parte es la regiÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â³n
                                                    region = parts[0]
                                                    value = parts[col_index + 1]
                                                    out_f.write(f'  YEAR_{year} {region} {value}\n')
                                        
                                        param_found = True
                                        break
                            
                            if not param_found:
                                logging.warning(f'Parameter {param_name} not found in year {year}')
                    
                    out_f.write(';\n\n')
                
            logging.info(f'Successfully processed misc data to {output_file}')
            
        except Exception as e:
            logging.error(f'Error processing misc data: {e}')
            import traceback
            traceback.print_exc()
        
    def _process_storage_dat(self):
        """
        Process storage power to energy data (storage_charge_time and storage_discharge_time)
        with year indexing
        """
        output_file = self.cs_dir / 'reg_storage_pathway.dat'
        
        try:
            with open(output_file, 'w') as out_f:
                out_f.write('# Storage power to energy data for pathway optimization\n')
                out_f.write('# Generated from year-specific reg_storage_power_to_energy.dat files\n\n')
                
                # Los parÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡metros a procesar
                storage_params = ['storage_charge_time', 'storage_discharge_time']
                
                # Procesar cada parÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡metro por separado
                for param_name in storage_params:
                    out_f.write(f'# {param_name} by year and region\n')
                    
                    # CAMBIO CRÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚ÂTICO: Eliminamos "default 0" de la declaraciÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â³n del parÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡metro
                    # Antes: out_f.write(f'param {param_name} default 0 :=\n')
                    # Ahora: Solo declaramos el parÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡metro sin especificar valor por defecto
                    out_f.write(f'param {param_name} :=\n')
                    
                    for year in self.years:
                        storage_file = self.data_dir / str(year) / 'reg_storage_power_to_energy.dat'
                        if storage_file.exists():
                            with open(storage_file, 'r') as f:
                                content = f.read()
                            
                            # Buscar el encabezado con mÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Âºltiples columnas
                            pattern = r'param\s*:\s*((?:storage_charge_time|storage_discharge_time)(?:\s+(?:storage_charge_time|storage_discharge_time))*)\s*:=(.*?);'
                            match = re.search(pattern, content, re.DOTALL)
                            
                            if match:
                                # Obtener los nombres de las columnas
                                param_names = match.group(1).split()
                                data_lines = match.group(2).strip().split('\n')
                                
                                # Encontrar el ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â­ndice de la columna para este parÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡metro
                                if param_name in param_names:
                                    col_index = param_names.index(param_name)
                                    
                                    # Procesar cada lÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â­nea de datos
                                    for line in data_lines:
                                        if line.strip() and not line.startswith('#'):
                                            parts = line.split()
                                            if len(parts) >= col_index + 3:  # RegiÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â³n + TecnologÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â­a + valores
                                                region = parts[0]
                                                technology = parts[1]
                                                value = parts[col_index + 2]  # +2 porque parts[0] es regiÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â³n, parts[1] es tecnologÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â­a
                                                
                                                # Escribir con formato: YEAR regiÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â³n tecnologÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â­a valor
                                                out_f.write(f'  YEAR_{year} {region} {technology} {value}\n')
                    
                    out_f.write(';\n\n')
                
            logging.info(f'Successfully processed storage data to {output_file}')
            
        except Exception as e:
            logging.error(f'Error processing storage data: {e}')
            import traceback
            traceback.print_exc()
        
    def _copy_pathway_files(self):
        """Copy pathway-specific .dat files from ALL_YEARS"""
        pathway_files = [
            'PES_seq_opti.dat',
            'PES_data_remaining.dat',
            'PES_data_all_years.dat',
            'PES_data_set_AGE_2021.dat',
            'PES_data_decom_allowed_2021.dat'
        ]
        
        for file_name in pathway_files:
            src_file = self.pathway_dir / file_name
            if src_file.exists():
                import shutil
                dest_file = self.cs_dir / file_name
                shutil.copy2(src_file, dest_file)
                logging.info(f'Copied {file_name} to output directory')
            else:
                logging.warning(f'Pathway file not found: {file_name}')
                
    def set_esom_pathway(self, ampl_path=None):
        """
        Set up the optimization model in AMPL
        
        This method configures AMPL with the model files and data files,
        following the same pattern as in esmc.py
        """
        # Model files - must be a list of Path objects
        mod_path = [
            self.project_dir / 'esmc' / 'energy_model' / 'ESMC_model_AMPL_pathway.mod',
            self.project_dir / 'esmc' / 'energy_model' / 'ESMC_obj_TotalCost_pathway.mod'
        ]
        
        # Data files - must be a list of all .dat file paths
        data_files = list(self.cs_dir.glob('*.dat'))
        data_path = [str(f) for f in data_files]  # Convert Path objects to strings
        
        logging.info(f'Setting up AMPL with {len(mod_path)} model files and {len(data_path)} data files')
        
        # Configure AMPL options (exactly as in esmc.py)
        cplex_options = [
            'baropt',
            'predual=-1',
            'barstart=4',
            'comptol=1e-4',
            'crossover=0',
            'timelimit 172800',
            'bardisplay=1',
            'display=2'
        ]
        
        cplex_options_str = ' '.join(cplex_options)
        
        ampl_options = {
            'show_stats': 3,
            'log_file': str(self.cs_dir / 'log.txt'),
            'presolve': 200,
            'times': 1,
            'gentimes': 1,
            'cplex_options': cplex_options_str
        }
        
        # Initialize optimizer using named arguments as in esmc.py
        self.esom = OptiProbl(
            mod_path=mod_path,           # List of .mod files
            data_path=data_path,         # List of .dat files
            options=ampl_options,        # Dictionary of options
            solver='cplex',              # Solver name
            ampl_path=ampl_path,         # Path to AMPL executable (optional)
            set_ampl=True                # Initialize AMPL immediately
        )
        
        logging.info('AMPL model configured successfully')
        
    def solve_esom_pathway(self):
        """Solve the optimization problem using AMPL"""
        logging.info('Starting pathway optimization...')
        # OptiProbl has a run_ampl() method that handles the solving
        self.esom.run_ampl()
        
    def get_results(self):
        """
        Extract key results from the optimization
        
        This method extracts the main results we're interested in:
        - F values (installed capacity by year, region, and technology)
        - Annual energy generation (calculated from F_t and t_op)
        
        The results are returned in a dictionary for easy access and processing.
        """
        results = {}
        
        # Get F values (installed capacity) using the get_var method from OptiProbl
        try:
            F_df = self.esom.get_var('F')
            results['F'] = F_df
            logging.info(f'Extracted F values: {len(F_df)} entries')
        except Exception as e:
            logging.error(f'Error extracting F: {e}')
            results['F'] = pd.DataFrame()
            
        # Get F_t values for calculating annual generation
        try:
            F_t_df = self.esom.get_var('F_t')
            results['F_t'] = F_t_df
            logging.info(f'Extracted F_t values: {len(F_t_df)} entries')
        except Exception as e:
            logging.error(f'Error extracting F_t: {e}')
            results['F_t'] = pd.DataFrame()
            
        # Get t_op values (operation time for each hour and typical day)
        try:
            # t_op is a parameter, not a variable, so we need to use get_param
            t_op_df = self.esom.get_param('t_op')
            results['t_op'] = t_op_df
            logging.info(f'Extracted t_op values: {len(t_op_df)} entries')
        except Exception as e:
            logging.error(f'Error extracting t_op: {e}')
            results['t_op'] = pd.DataFrame()
            
        # Calculate annual generation if we have both F_t and t_op
        if not results['F_t'].empty and not results['t_op'].empty:
            try:
                # Annual generation = sum over all hours and typical days of (F_t * t_op)
                # This calculation depends on the exact structure of your data
                # You may need to adjust based on how the indices are structured
                logging.info('Calculating annual energy generation...')
                
                # Note: The exact calculation will depend on your model structure
                # This is a placeholder that you'll need to adapt
                annual_gen = self._calculate_annual_generation(results['F_t'], results['t_op'])
                results['annual_generation'] = annual_gen
                
            except Exception as e:
                logging.error(f'Error calculating annual generation: {e}')
                results['annual_generation'] = pd.DataFrame()
        
        return results
    
    def _calculate_annual_generation(self, F_t_df, t_op_df):
        """
        Calculate annual energy generation from F_t and t_op
        
        This is a placeholder method that you'll need to implement based on
        your specific model structure and how the data is indexed.
        
        The formula should implement:
        sum {t in PERIODS, h in HOUR_OF_PERIOD[t], td in TYPICAL_DAY_OF_PERIOD[t]} 
            (F_t[y,c,j,h,td] * t_op[h,td])
        """
        # TODO: Implement the actual calculation based on your model structure
        logging.warning('Annual generation calculation not yet implemented - returning empty DataFrame')
        return pd.DataFrame()    
    def _td_to_year_with_top(self, df_td, value_col):
        """
        Convert typical day values to yearly values using t_op weights
        
        FIXED: Now properly expands typical days to all 8760 hours before summing.
        
        This method replicates the AMPL calculation:
            sum {t in PERIODS, h in HOUR_OF_PERIOD[t], td in TYPICAL_DAY_OF_PERIOD[t]} 
                value[h,td] * t_op[h,td]
        
        Parameters:
        -----------
        df_td : DataFrame
            DataFrame with MultiIndex including 'Typical_days' and 'Hours' levels
            Shape: (n_TDs * 24 * n_entities) rows
        value_col : str
            Name of the column containing values to aggregate
            
        Returns:
        --------
        DataFrame with yearly aggregated values
            Shape: (n_entities) rows with annual sums
        """
        try:
            # Step 1: Ensure t_h_td mapping is available
            if self.t_h_td is None:
                logging.info("Generating t_h_td mapping for temporal aggregation")
                self.generate_t_h_td()
            
            if self.t_h_td is None:
                raise ValueError("t_h_td mapping could not be generated")
            
            # Step 2: Reset index to work with columns
            df_reset = df_td.reset_index()
            
            # Step 3: Identify which columns are TD and Hour (they may have different names)
            td_col = None
            h_col = None
            
            for col in df_reset.columns:
                col_lower = str(col).lower()
                if 'typical' in col_lower or col == 'Typical_days':
                    td_col = col
                elif 'hour' in col_lower and 'of_y' not in col_lower and col != 'Hours':
                    h_col = col
                elif col == 'Hours':
                    h_col = col
            
            if td_col is None or h_col is None:
                logging.error(f"Could not identify TD and Hour columns. Available: {df_reset.columns.tolist()}")
                raise ValueError("Could not identify Typical_days and Hours columns")
            
            # Step 4: Prepare data for merging - rename TD column to match t_h_td
            df_for_merge = df_reset.copy()
            
            # Map to TD_number if needed
            if td_col != 'TD_number':
                unique_tds = sorted(df_for_merge[td_col].unique())
                td_mapping = {td: idx+1 for idx, td in enumerate(unique_tds)}
                df_for_merge['TD_number'] = df_for_merge[td_col].map(td_mapping)
            else:
                df_for_merge['TD_number'] = df_for_merge[td_col]
            
            # Rename hour column to H_of_D for merging
            if h_col != 'H_of_D':
                df_for_merge['H_of_D'] = df_for_merge[h_col]
            
            # Step 5: CRITICAL - Expand from typical days to all 8760 hours
            logging.debug(f"Expanding from {len(df_for_merge)} typical day hours to 8760 year hours")
            
            group_cols = [col for col in df_for_merge.columns 
                         if col not in [td_col, h_col, 'TD_number', 'H_of_D', value_col]]
            
            # Merge with t_h_td to expand to full year
            df_expanded = self.t_h_td[['TD_number', 'H_of_D', 'H_of_Y']].merge(
                df_for_merge,
                on=['TD_number', 'H_of_D'],
                how='left'
            )
            
            if len(df_expanded) == 0:
                raise ValueError("Merge with t_h_td produced empty result")
            
            # Step 6: Get t_op parameter from AMPL
            t_op = self.esom.get_param('t_op')
            t_op_reset = t_op.reset_index()
            
            # Normalize t_op column names
            rename_map = {}
            for col in t_op_reset.columns:
                if 'HOUR' in col.upper() and col not in ['Hours', 'H_of_D']:
                    rename_map[col] = 'Hours'
                elif 'TYPICAL' in col.upper() and 'DAY' in col.upper() and col not in ['Typical_days', 'TD_number']:
                    rename_map[col] = 'Typical_days'
            
            if rename_map:
                t_op_reset.rename(columns=rename_map, inplace=True)
            
            # Create TD_number mapping for t_op if needed
            if 'TD_number' not in t_op_reset.columns:
                if 'Typical_days' in t_op_reset.columns:
                    unique_tds_top = sorted(t_op_reset['Typical_days'].unique())
                    td_mapping_top = {td: idx+1 for idx, td in enumerate(unique_tds_top)}
                    t_op_reset['TD_number'] = t_op_reset['Typical_days'].map(td_mapping_top)
            
            if 'H_of_D' not in t_op_reset.columns:
                if 'Hours' in t_op_reset.columns:
                    t_op_reset['H_of_D'] = t_op_reset['Hours']
            
            # Step 7: Merge t_op with expanded data
            df_with_weights = df_expanded.merge(
                t_op_reset[['TD_number', 'H_of_D', 't_op']],
                on=['TD_number', 'H_of_D'],
                how='left'
            )
            
            # Handle missing t_op values
            if df_with_weights['t_op'].isna().all():
                logging.warning("t_op merge failed, using t_op=1 for all hours")
                df_with_weights['t_op'] = 1.0
            elif df_with_weights['t_op'].isna().any():
                n_missing = df_with_weights['t_op'].isna().sum()
                logging.warning(f"{n_missing} t_op values are NaN, filling with 1.0")
                df_with_weights['t_op'].fillna(1.0, inplace=True)
            
            # Step 8: Calculate weighted values
            df_with_weights[f'{value_col}_weighted'] = df_with_weights[value_col] * df_with_weights['t_op']
            
            # Step 9: Sum over all 8760 hours
            if group_cols:
                result = df_with_weights.groupby(group_cols)[f'{value_col}_weighted'].sum().to_frame()
                result.columns = [value_col]
            else:
                total = df_with_weights[f'{value_col}_weighted'].sum()
                result = pd.DataFrame({value_col: [total]})
            
            # Log statistics
            total_expanded = df_with_weights[f'{value_col}_weighted'].sum()
            total_td_only = df_for_merge[value_col].sum()
            scaling = total_expanded/total_td_only if total_td_only > 0 else 0
            
            logging.info(f"TD to year conversion: TD sum={total_td_only:.2f}, "
                        f"Annual sum={total_expanded:.2f}, Scaling={scaling:.2f}x")
            
            return result
            
        except Exception as e:
            logging.error(f"Error in _td_to_year_with_top: {e}")
            import traceback
            traceback.print_exc()
            
            # Fallback: try to use from_td_to_year
            try:
                logging.warning("Attempting fallback to from_td_to_year method")
                df_reset = df_td.reset_index()
                
                # Find TD and Hour columns
                td_col = None
                h_col = None
                for col in df_reset.columns:
                    col_lower = str(col).lower()
                    if 'typical' in col_lower:
                        td_col = col
                    elif 'hour' in col_lower and 'of_y' not in col_lower:
                        h_col = col
                
                if td_col and h_col:
                    df_for_conversion = df_reset.set_index([td_col, h_col])
                    expanded = self.from_td_to_year(df_for_conversion)
                    
                    # Apply t_op weighting
                    t_op = self.esom.get_param('t_op')
                    # For simplicity, just sum without detailed t_op weighting
                    group_cols = [col for col in expanded.columns if col not in [value_col]]
                    if group_cols:
                        result = expanded.groupby(group_cols)[value_col].sum().to_frame()
                    else:
                        result = pd.DataFrame({value_col: [expanded[value_col].sum()]})
                    
                    logging.info("Fallback successful using from_td_to_year")
                    return result
            except Exception as e2:
                logging.error(f"Fallback also failed: {e2}")
            
            # Final fallback: direct sum (INCORRECT)
            logging.error("All methods failed, using direct sum (WILL BE ~22x TOO SMALL!)")
            group_cols = [name for name in df_td.index.names if name not in ['Hours', 'Typical_days']]
            if group_cols:
                return df_td.groupby(group_cols)[value_col].sum().to_frame()
            else:
                return pd.DataFrame({value_col: [df_td[value_col].sum()]})

    def solve_esom_pathway(self):
        """Solve the optimization problem using AMPL"""
        logging.info('Starting pathway optimization...')
        self.esom.run_ampl()
    
    def get_total_cost_pathway(self):
        """Get the total annualized cost for different regions and years"""
        logging.info('Getting TotalCost (pathway)')
        
        total_cost = self.esom.get_var('TotalCost').reset_index()
        
        if 'YEARS' in total_cost.columns:
            total_cost.rename(columns={'YEARS': 'Years'}, inplace=True)
        if 'REGIONS' in total_cost.columns:
            total_cost.rename(columns={'REGIONS': 'Regions'}, inplace=True)
        
        if hasattr(self, 'regions_names'):
            regions_names = self.regions_names
        else:
            regions_names = sorted(total_cost['Regions'].unique())
        
        total_cost['Years'] = pd.Categorical(total_cost['Years'], 
                                             sorted(total_cost['Years'].unique()))
        total_cost['Regions'] = pd.Categorical(total_cost['Regions'], regions_names)
        
        total_cost = total_cost.set_index(['Years', 'Regions'])
        total_cost.sort_index(inplace=True)
        
        self.results['TotalCost'] = total_cost
        self.results_all['TotalCost'] = total_cost.groupby('Years', observed=True).sum()
        
        return
    
    def get_cost_breakdown_pathway(self):
        """Gets cost breakdown with YEARS dimension
        
        DEFINITIVO: Extracts tau directly from AMPL (already calculated as parameter)
        This matches the structure of both one-year and pathway AMPL models where:
        - C_inv = c_inv * F (without annualization)
        - tau is defined in AMPL as parameter
        - Cost breakdown needs C_inv * tau for annualized investment cost
        """
        logging.info('Getting Cost_breakdown (pathway)')
        
        # Get the different cost variables
        c_inv = self.esom.get_var('C_inv')
        c_maint = self.esom.get_var('C_maint')
        c_op = self.esom.get_var('C_op')
        
        # Rename indices
        for df in [c_inv, c_maint, c_op]:
            if 'YEARS' in df.index.names:
                df.index = df.index.set_names('Years', level='YEARS')
            if 'REGIONS' in df.index.names:
                df.index = df.index.set_names('Regions', level='REGIONS')
            # For technologies/resources
            for level_name in df.index.names:
                if 'TECHNOLOGIES' in str(level_name).upper():
                    df.index = df.index.set_names('Elements', level=level_name)
                elif 'RESOURCES' in str(level_name).upper():
                    df.index = df.index.set_names('Elements', level=level_name)
        
        # Get years list
        years_list = sorted(c_inv.index.get_level_values('Years').unique())
        
        # Get regions list
        if hasattr(self, 'regions_names'):
            regions_names = self.regions_names
        else:
            regions_names = sorted(c_inv.index.get_level_values('Regions').unique())
        
        # EXTRACT TAU DIRECTLY FROM AMPL 
        # In AMPL: param tau {y in YEARS, c in REGIONS, i in TECHNOLOGIES} := 
        #              i_rate * (1 + i_rate)^lifetime[y,c,i] / (((1 + i_rate)^lifetime[y,c,i]) - 1);
        try:
            tau_param = self.esom.get_param('tau')
            
            # Rename indices to match c_inv
            if 'YEARS' in tau_param.index.names:
                tau_param.index = tau_param.index.set_names('Years', level='YEARS')
            if 'REGIONS' in tau_param.index.names:
                tau_param.index = tau_param.index.set_names('Regions', level='REGIONS')
            for level_name in tau_param.index.names:
                if 'TECHNOLOGIES' in str(level_name).upper():
                    tau_param.index = tau_param.index.set_names('Elements', level=level_name)
            
            # Annualize C_inv
            if isinstance(tau_param, pd.DataFrame):
                # Find the column with tau values (should be named 'tau' or be the only numeric column)
                value_col = None
                if 'tau' in tau_param.columns:
                    value_col = 'tau'
                else:
                    # Find the first numeric column
                    numeric_cols = tau_param.select_dtypes(include=[np.number]).columns
                    if len(numeric_cols) > 0:
                        value_col = numeric_cols[0]
                
                if value_col is not None:
                    # Convert to Series with proper index
                    tau_param = tau_param[value_col]
                else:
                    # If we can't find the value column, try to use the whole DataFrame
                    # but drop any columns that look like index columns
                    cols_to_drop = []
                    for col in tau_param.columns:
                        col_lower = str(col).lower()
                        if any(x in col_lower for x in ['year', 'region', 'technolog', 'element', ' in ']):
                            cols_to_drop.append(col)
                    if cols_to_drop:
                        logging.warning(f"Dropping extra columns from tau_param: {cols_to_drop}")
                        tau_param = tau_param.drop(columns=cols_to_drop)
                    # If only one column left, convert to Series
                    if len(tau_param.columns) == 1:
                        tau_param = tau_param.iloc[:, 0]
            
            # Now tau_param should be a clean Series indexed by (Years, Regions, Elements)
            # Multiply c_inv by tau element-wise (same as one-year model approach)
            c_inv_ann = c_inv.mul(tau_param, axis=0)
            
            logging.info("Successfully extracted tau parameter from AMPL and annualized C_inv")
            
        except Exception as e:
            logging.error(f"CRITICAL: Could not extract tau from AMPL: {e}")
            logging.error("C_inv will be incorrect (not annualized). Check AMPL model has tau parameter.")
            # Use C_inv as-is (will be incorrect)
            c_inv_ann = c_inv
        
        # Merge costs into cost breakdown
        cost_breakdown = c_inv_ann.merge(c_maint, left_index=True, right_index=True, how='outer') \
            .merge(c_op, left_index=True, right_index=True, how='outer')
        
        # Set categorical data for sorting
        cost_breakdown = cost_breakdown.reset_index()
        cost_breakdown['Years'] = pd.Categorical(cost_breakdown['Years'], years_list)
        cost_breakdown['Regions'] = pd.Categorical(cost_breakdown['Regions'], regions_names)
        
        # Get elements order
        elements_order = sorted(cost_breakdown['Elements'].unique())
        cost_breakdown['Elements'] = pd.Categorical(cost_breakdown['Elements'], elements_order)
        
        cost_breakdown.sort_values(by=['Years', 'Regions', 'Elements'], axis=0, ignore_index=True, inplace=True)
        cost_breakdown.set_index(['Years', 'Regions', 'Elements'], inplace=True)
        
        # Put very small values as nan
        threshold = 1e-2
        # Select only numeric columns for masking
        numeric_cols = cost_breakdown.select_dtypes(include=[np.number]).columns
        if len(numeric_cols) > 0:
            cost_breakdown[numeric_cols] = cost_breakdown[numeric_cols].mask(
                (cost_breakdown[numeric_cols] > -threshold) & (cost_breakdown[numeric_cols] < threshold), 
                np.nan
            )
        
        # Store into results
        self.results['Cost_breakdown'] = cost_breakdown
        self.results_all['Cost_breakdown'] = cost_breakdown.groupby(['Years', 'Elements'], observed=True).sum()
        
        return
    
    def get_gwp_breakdown_pathway(self):
        """Gets GWP breakdown with YEARS dimension
        
        DEFINITIVO: Properly uses lifetime from AMPL for annualization
        """
        logging.info('Getting Gwp_breakdown (pathway)')
        
        # Get GWP variables
        gwp_constr = self.esom.get_var('GWP_constr')
        gwp_op = self.esom.get_var('GWP_op')
        
        try:
            co2_net = self.esom.get_var('CO2_net')
        except:
            co2_net = None
        
        # Rename indices
        for df in [gwp_constr, gwp_op]:
            if df is None:
                continue
            if 'YEARS' in df.index.names:
                df.index = df.index.set_names('Years', level='YEARS')
            if 'REGIONS' in df.index.names:
                df.index = df.index.set_names('Regions', level='REGIONS')
            for level_name in df.index.names:
                if 'TECHNOLOGIES' in str(level_name).upper():
                    df.index = df.index.set_names('Elements', level=level_name)
                elif 'RESOURCES' in str(level_name).upper():
                    df.index = df.index.set_names('Elements', level=level_name)
        
        if co2_net is not None:
            if 'YEARS' in co2_net.index.names:
                co2_net.index = co2_net.index.set_names('Years', level='YEARS')
            if 'REGIONS' in co2_net.index.names:
                co2_net.index = co2_net.index.set_names('Regions', level='REGIONS')
            for level_name in co2_net.index.names:
                if 'RESOURCES' in str(level_name).upper():
                    co2_net.index = co2_net.index.set_names('Elements', level=level_name)
        
        # Get years and regions
        years_list = sorted(gwp_constr.index.get_level_values('Years').unique())
        if hasattr(self, 'regions_names'):
            regions_names = self.regions_names
        else:
            regions_names = sorted(gwp_constr.index.get_level_values('Regions').unique())
        
        # Get lifetime to annualize GWP_constr
        try:
            lifetime_param = self.esom.get_param('lifetime')
            
            # Rename indices
            if 'YEARS' in lifetime_param.index.names:
                lifetime_param.index = lifetime_param.index.set_names('Years', level='YEARS')
            if 'REGIONS' in lifetime_param.index.names:
                lifetime_param.index = lifetime_param.index.set_names('Regions', level='REGIONS')
            for level_name in lifetime_param.index.names:
                if 'TECHNOLOGIES' in str(level_name).upper():
                    lifetime_param.index = lifetime_param.index.set_names('Elements', level=level_name)
            
            # Annualize GWP_constr by dividing by lifetime
            gwp_constr_ann = gwp_constr.merge(lifetime_param, left_index=True, right_index=True, how='left')
            gwp_constr_ann['GWP_constr'] = gwp_constr_ann['GWP_constr'] / gwp_constr_ann['lifetime']
            gwp_constr_ann = gwp_constr_ann[['GWP_constr']]
            
            logging.info("Successfully used lifetime from AMPL for GWP annualization")
            
        except Exception as e:
            logging.error(f"Could not get lifetime from AMPL for GWP: {e}. GWP_constr will not be annualized.")
            gwp_constr_ann = gwp_constr
        
        # Merge emissions into gwp_breakdown
        if co2_net is not None:
            gwp_breakdown = gwp_constr_ann.merge(gwp_op, left_index=True, right_index=True, how='outer') \
                .merge(co2_net, left_index=True, right_index=True, how='outer')
        else:
            gwp_breakdown = gwp_constr_ann.merge(gwp_op, left_index=True, right_index=True, how='outer')
        
        # Set categorical data for sorting
        gwp_breakdown = gwp_breakdown.reset_index()
        gwp_breakdown['Years'] = pd.Categorical(gwp_breakdown['Years'], years_list)
        gwp_breakdown['Regions'] = pd.Categorical(gwp_breakdown['Regions'], regions_names)
        
        elements_order = sorted(gwp_breakdown['Elements'].unique())
        gwp_breakdown['Elements'] = pd.Categorical(gwp_breakdown['Elements'], elements_order)
        
        gwp_breakdown.sort_values(by=['Years', 'Regions', 'Elements'], axis=0, ignore_index=True, inplace=True)
        gwp_breakdown.set_index(['Years', 'Regions', 'Elements'], inplace=True)
        
        # Put very small values as nan
        threshold = 1e-5
        # Select only numeric columns for masking
        numeric_cols = gwp_breakdown.select_dtypes(include=[np.number]).columns
        if len(numeric_cols) > 0:
            gwp_breakdown[numeric_cols] = gwp_breakdown[numeric_cols].mask(
                (gwp_breakdown[numeric_cols] > -threshold) & (gwp_breakdown[numeric_cols] < threshold), 
                np.nan
            )
        
        # Store into results
        self.results['Gwp_breakdown'] = gwp_breakdown
        self.results_all['Gwp_breakdown'] = gwp_breakdown.groupby(['Years', 'Elements'], observed=True).sum()
        
        return
    
    def get_resources_and_exchanges_pathway(self, save_hourly:list=[]):
        """Gets yearly use of resources and exchanges with YEARS dimension
        
        DEFINITIVO: Fixed merge issue with explicit DataFrame creation
        """
        logging.info('Getting Resources and Exchanges (pathway)')
        
        # Get R_t_local and R_t_exterior
        r_t_local = self.esom.get_var('R_t_local')
        r_t_exterior = self.esom.get_var('R_t_exterior')
        
        # Rename indices properly
        for df in [r_t_local, r_t_exterior]:
            if 'YEARS' in df.index.names:
                df.index = df.index.set_names('Years', level='YEARS')
            if 'REGIONS' in df.index.names:
                df.index = df.index.set_names('Regions', level='REGIONS')
            if 'I in RESOURCES' in df.index.names:
                df.index = df.index.set_names('Resources', level='I in RESOURCES')
            if 'HOURS' in df.index.names:
                df.index = df.index.set_names('Hours', level='HOURS')
            if 'TYPICAL_DAYS' in df.index.names:
                df.index = df.index.set_names('Typical_days', level='TYPICAL_DAYS')
        
        # Get region and year lists
        if hasattr(self, 'regions_names'):
            regions_names = self.regions_names
        else:
            regions_names = sorted(r_t_local.index.get_level_values('Regions').unique())
        
        years_list = sorted(r_t_local.index.get_level_values('Years').unique())
        
        # Calculate yearly values using t_op weighting
        resources_local_list = []
        resources_exterior_list = []
        
        for year in years_list:
            # Extract data for this year
            r_t_local_year = r_t_local.xs(year, level='Years')
            r_t_exterior_year = r_t_exterior.xs(year, level='Years')
            
            # Use from_td_to_year to expand to 8760 hours, then sum
            if self.ta is not None:
                r_year_local = self.ta.from_td_to_year(
                    ts_td=r_t_local_year.reset_index().set_index(['Typical_days', 'Hours'])
                ).groupby(['Regions', 'Resources']).sum().rename(columns={'R_t_local': 'R_year_local'})
                
                r_year_exterior = self.ta.from_td_to_year(
                    ts_td=r_t_exterior_year.reset_index().set_index(['Typical_days', 'Hours'])
                ).groupby(['Regions', 'Resources']).sum().rename(columns={'R_t_exterior': 'R_year_exterior'})
            else:
                logging.warning("No temporal aggregation available, using direct TD sum")
                r_year_local = r_t_local_year.groupby(['Regions', 'Resources']).sum().rename(columns={'R_t_local': 'R_year_local'})
                r_year_exterior = r_t_exterior_year.groupby(['Regions', 'Resources']).sum().rename(columns={'R_t_exterior': 'R_year_exterior'})

            # Add year dimension back
            r_year_local = r_year_local.reset_index()
            r_year_local['Years'] = year
            r_year_local = r_year_local.set_index(['Years', 'Regions', 'Resources'])

            r_year_exterior = r_year_exterior.reset_index()
            r_year_exterior['Years'] = year
            r_year_exterior = r_year_exterior.set_index(['Years', 'Regions', 'Resources'])
            r_year_exterior.columns = ['R_year_exterior']

            resources_local_list.append(r_year_local)
            resources_exterior_list.append(r_year_exterior)
        
        # Combine all years
        r_year_local_all = pd.concat(resources_local_list)
        r_year_exterior_all = pd.concat(resources_exterior_list)
        
        # Merge local and exterior
        resources = r_year_local_all.merge(r_year_exterior_all, left_index=True, right_index=True, how='outer')
        resources = resources.fillna(0)
        
        # ===== EXCHANGES AND TRANSFER CAPACITY CALCULATIONS =====
        # Get list of resources exchanged
        try:
            network_exch_r = self.esom.ampl.get_set('EXCHANGE_NETWORK_R').getValues().toList()
        except:
            network_exch_r = []
        
        try:
            freight_exch_r = self.esom.ampl.get_set('EXCHANGE_FREIGHT_R').getValues().toList()
        except:
            freight_exch_r = []
        
        r_exch = network_exch_r.copy()
        r_exch.extend(freight_exch_r)
        
        if len(r_exch) > 0:
            # Get exchange variables
            try:
                exch_imp = self.esom.get_var('Exch_imp')
                exch_exp = self.esom.get_var('Exch_exp')
                
                # Rename indices
                for df in [exch_imp, exch_exp]:
                    if 'YEARS' in df.index.names:
                        df.index = df.index.set_names('Years', level='YEARS')
                    # The index structure for Exch_imp and Exch_exp is: YEARS, From, To, Resources, Hours, Typical_days
                    index_names = list(df.index.names)
                    new_names = []
                    regions_count = 0
                    for name in index_names:
                        if 'REGIONS' in str(name).upper():
                            if regions_count == 0:
                                new_names.append('From')
                                regions_count += 1
                            else:
                                new_names.append('To')
                        elif 'EXCHANGE_R' in str(name).upper() or ('RESOURCES' in str(name).upper() and 'From' in new_names):
                            new_names.append('Resources')
                        elif 'HOURS' in str(name).upper():
                            new_names.append('Hours')
                        elif 'TYPICAL' in str(name).upper():
                            new_names.append('Typical_days')
                        else:
                            new_names.append(name)
                    df.index.names = new_names
                
                # Get Transfer_capacity
                transfer_capacity = self.esom.get_var('Transfer_capacity')
                if 'YEARS' in transfer_capacity.index.names:
                    transfer_capacity.index = transfer_capacity.index.set_names('Years', level='YEARS')
                # Rename other indices
                index_names = list(transfer_capacity.index.names)
                new_names = []
                regions_count = 0
                for name in index_names:
                    if 'REGIONS' in str(name).upper():
                        if regions_count == 0:
                            new_names.append('From')
                            regions_count += 1
                        else:
                            new_names.append('To')
                    elif 'EXCHANGE_NETWORK_R' in str(name).upper() or ('RESOURCES' in str(name).upper() and 'From' in new_names):
                        new_names.append('Resources')
                    elif 'NETWORK_TYPE' in str(name).upper():
                        new_names.append('Network_type')
                    else:
                        new_names.append(name)
                transfer_capacity.index.names = new_names
                
                # Store all transfer capacity
                all_tc = transfer_capacity.copy()
                threshold = 1e-3
                # Select only numeric columns for masking
                numeric_cols = all_tc.select_dtypes(include=[np.number]).columns
                if len(numeric_cols) > 0:
                    all_tc[numeric_cols] = all_tc[numeric_cols].mask(
                        (all_tc[numeric_cols] > -threshold) & (all_tc[numeric_cols] < threshold), 
                        np.nan
                    )
                
                # Process exchanges for each year
                exchanges_year_list = []
                exchanges_year_all_list = []
                r_year_import_list = []
                r_year_export_list = []
                
                for year in years_list:
                    # Extract exchanges for this year
                    exch_imp_year = exch_imp.xs(year, level='Years')
                    exch_exp_year = exch_exp.xs(year, level='Years')
                    
                    # FIXED: Clean exchanges from double fictive fluxes without creating duplicate columns
                    # Convert both to DataFrame with explicit column names
                    exch_imp_df = pd.DataFrame({'Exch_imp': exch_imp_year['Exch_imp']})
                    exch_exp_df = pd.DataFrame({'Exch_exp': exch_exp_year['Exch_exp']})
                    
                    # Merge them
                    exch = exch_exp_df.merge(exch_imp_df, left_index=True, right_index=True, how='outer')
                    exch['Balance'] = exch['Exch_exp'] - exch['Exch_imp']
                    
                    # Replace Exch_imp and Exch_exp by values deduced from Balance
                    threshold = 1e-6
                    exch['Exch_imp'] = exch['Balance'].mask((exch['Balance'] > -threshold), np.nan)
                    exch['Exch_exp'] = exch['Balance'].mask((exch['Balance'] < threshold), np.nan)
                    exch['Balance'] = exch['Balance'].mask((exch['Balance'].abs() < threshold), np.nan)
                    
                    # Compute total over the year for each LINK
                    exch_for_agg = exch[['Exch_exp', 'Exch_imp']].copy()

                    if self.ta is not None:
                        # Expand both columns at once
                        exch_expanded = self.ta.from_td_to_year(
                            ts_td=exch_for_agg.reset_index().set_index(['Typical_days', 'Hours'])
                        ).groupby(['From', 'To', 'Resources']).sum()
                        
                        exchanges_year_data = exch_expanded
                    else:
                        logging.warning("No temporal aggregation available, using direct TD sum")
                        exchanges_year_data = exch_for_agg.groupby(['From', 'To', 'Resources']).sum()
                    
                    # Compute total over the year for each REGION
                    r_exch_region = exchanges_year_data.reset_index().groupby(['From', 'Resources']).sum(numeric_only=True).abs()
                    
                    # Keep only one direction per link (Exch_exp)
                    exchanges_year_data = exchanges_year_data[['Exch_exp']].rename(columns={'Exch_exp': 'Exchanges_year'})
                    
                    # Compute utilization factor of lines
                    transfer_capacity_year = transfer_capacity.xs(year, level='Years')
                    transfer_capacity_grouped = transfer_capacity_year.groupby(['From', 'To', 'Resources']).sum()
                    
                    exchanges_year_data = exchanges_year_data.merge(
                        transfer_capacity_grouped, 
                        how='outer', 
                        left_index=True, 
                        right_index=True
                    )
                    exchanges_year_data['Utilization_factor'] = \
                        exchanges_year_data['Exchanges_year'] / (transfer_capacity_grouped['Transfer_capacity'] * 8760)
                    
                    # Add year dimension back
                    exchanges_year_data = exchanges_year_data.assign(Years=year).reset_index()
                    exchanges_year_data = exchanges_year_data.set_index(['Years', 'From', 'To', 'Resources'])
                    
                    exchanges_year_list.append(exchanges_year_data)
                    
                    # Compute balance by REGION at each hour for total fluxes
                    r_t = exch.groupby(['Resources', 'From', 'Hours', 'Typical_days']).sum()
                    # Keep only net exporter
                    r_t['Balance'] = r_t['Balance'].mask(r_t['Balance'] < 0, 0)
                    # Sum over all regions and all hours to get total fluxes
                    r_t = r_t.groupby(['Resources', 'Hours', 'Typical_days']).sum()
                    
                    if self.ta is not None:
                        exchanges_year_all_data = self.ta.from_td_to_year(
                            ts_td=r_t[['Balance']].reset_index().set_index(['Typical_days', 'Hours'])
                        ).groupby(['Resources']).sum().rename(columns={'Balance': 'Exchanges_year'})
                    else:
                        logging.warning("No temporal aggregation available, using direct TD sum")
                        exchanges_year_all_data = r_t[['Balance']].groupby(['Resources']).sum().rename(columns={'Balance': 'Exchanges_year'})
                    exchanges_year_all_data = exchanges_year_all_data.assign(Years=year).reset_index()
                    exchanges_year_all_data = exchanges_year_all_data.set_index(['Years', 'Resources'])
                    exchanges_year_all_list.append(exchanges_year_all_data)
                    
                    # Compute R_year_import and R_year_export from exchanges
                    r_year_export = pd.DataFrame(r_exch_region['Exch_exp']).rename(columns={'Exch_exp': 'R_year_export'})
                    r_year_export.index.names = ['Regions', 'Resources']
                    r_year_export = r_year_export.assign(Years=year).reset_index()
                    r_year_export = r_year_export.set_index(['Years', 'Regions', 'Resources'])
                    r_year_export_list.append(r_year_export)
                    
                    r_year_import = pd.DataFrame(r_exch_region['Exch_imp']).rename(columns={'Exch_imp': 'R_year_import'})
                    r_year_import.index.names = ['Regions', 'Resources']
                    r_year_import = r_year_import.assign(Years=year).reset_index()
                    r_year_import = r_year_import.set_index(['Years', 'Regions', 'Resources'])
                    r_year_import_list.append(r_year_import)
                
                # Combine all years
                if exchanges_year_list:
                    exchanges_year = pd.concat(exchanges_year_list)
                    
                    # Set categorical data for sorting
                    exchanges_year = exchanges_year.reset_index()
                    exchanges_year['Years'] = pd.Categorical(exchanges_year['Years'], years_list)
                    exchanges_year['From'] = pd.Categorical(exchanges_year['From'], regions_names)
                    exchanges_year['To'] = pd.Categorical(exchanges_year['To'], regions_names)
                    resource_order = sorted(exchanges_year['Resources'].unique())
                    exchanges_year['Resources'] = pd.Categorical(exchanges_year['Resources'], resource_order)
                    exchanges_year.sort_values(by=['Years', 'From', 'To', 'Resources'], axis=0, ignore_index=True, inplace=True)
                    exchanges_year.set_index(['Years', 'From', 'To', 'Resources'], inplace=True)
                    
                    # Put very small values as nan
                    threshold = 1e-3
                    # Select only numeric columns for masking
                    numeric_cols = exchanges_year.select_dtypes(include=[np.number]).columns
                    if len(numeric_cols) > 0:
                        exchanges_year[numeric_cols] = exchanges_year[numeric_cols].mask(
                            (exchanges_year[numeric_cols] > -threshold) & (exchanges_year[numeric_cols] < threshold), 
                            np.nan
                        )
                else:
                    exchanges_year = pd.DataFrame()
                
                if exchanges_year_all_list:
                    exchanges_year_all = pd.concat(exchanges_year_all_list)
                else:
                    exchanges_year_all = pd.DataFrame()
                
                # Merge import/export into resources
                if r_year_import_list:
                    r_year_import_all = pd.concat(r_year_import_list)
                    r_year_export_all = pd.concat(r_year_export_list)
                    
                    resources = resources.merge(r_year_import_all, left_index=True, right_index=True, how='outer')
                    resources = resources.merge(r_year_export_all, left_index=True, right_index=True, how='outer')
                
                logging.info(f"Successfully calculated exchanges and transfer capacity")
                
            except Exception as e:
                logging.error(f"Could not calculate exchanges and transfer capacity: {e}")
                import traceback
                traceback.print_exc()
                exchanges_year = pd.DataFrame()
                exchanges_year_all = pd.DataFrame()
                all_tc = pd.DataFrame()
        else:
            # No exchanges
            exchanges_year = pd.DataFrame()
            exchanges_year_all = pd.DataFrame()
            all_tc = pd.DataFrame()
        
        # Set categorical data for sorting resources
        resources = resources.reset_index()
        resources['Years'] = pd.Categorical(resources['Years'], years_list)
        resources['Regions'] = pd.Categorical(resources['Regions'], regions_names)
        
        resource_order = sorted(resources['Resources'].unique())
        resources['Resources'] = pd.Categorical(resources['Resources'], resource_order)
        
        resources.sort_values(by=['Years', 'Regions', 'Resources'], ignore_index=True, inplace=True)
        resources.set_index(['Years', 'Regions', 'Resources'], inplace=True)
        
        # Put very small values as nan
        threshold = 1e-5
        # Select only numeric columns for masking
        numeric_cols = resources.select_dtypes(include=[np.number]).columns
        if len(numeric_cols) > 0:
            resources[numeric_cols] = resources[numeric_cols].mask(
                (resources[numeric_cols] > -threshold) & (resources[numeric_cols] < threshold), 
                np.nan
            )
        
        # Store results
        self.results['Resources'] = resources
        self.results_all['Resources'] = resources.groupby(['Years', 'Resources'], observed=True).sum()
        
        self.results['Exchanges_year'] = exchanges_year
        if not exchanges_year.empty:
            self.results_all['Exchanges_year'] = exchanges_year_all
        else:
            self.results_all['Exchanges_year'] = pd.DataFrame()
            
        self.results['Transfer_capacity'] = all_tc
        if not all_tc.empty:
            # For results_all, group by Years, Resources, Network_type and divide by 2 (bidirectional)
            all_tc_reset = all_tc.reset_index()
            if 'Network_type' in all_tc_reset.columns:
                self.results_all['Transfer_capacity'] = all_tc_reset.groupby(['Years', 'Resources', 'Network_type'])['Transfer_capacity'].sum().to_frame() / 2
            else:
                self.results_all['Transfer_capacity'] = all_tc_reset.groupby(['Years', 'Resources'])['Transfer_capacity'].sum().to_frame() / 2
        else:
            self.results_all['Transfer_capacity'] = pd.DataFrame()
        
        # Save hourly if requested
        if 'Resources' in save_hourly:
            try:
                # Combine local and exterior for hourly data
                r_t_local_hourly = r_t_local.reset_index().set_index(
                    ['Years', 'Regions', 'Resources', 'Typical_days', 'Hours']).sort_index()
                r_t_exterior_hourly = r_t_exterior.reset_index().set_index(
                    ['Years', 'Regions', 'Resources', 'Typical_days', 'Hours']).sort_index()
                
                hourly_resources = r_t_local_hourly.merge(r_t_exterior_hourly, left_index=True, right_index=True, how='outer')
                self.hourly_results['Resources'] = hourly_resources
            except Exception as e:
                logging.warning(f"Could not save hourly Resources: {e}")
        
        if 'Exchanges' in save_hourly and len(r_exch) > 0:
            try:
                # Save exchanges hourly data
                exch_all_years = []
                for year in years_list:
                    exch_imp_year = exch_imp.xs(year, level='Years')
                    exch_exp_year = exch_exp.xs(year, level='Years')
                    
                    exch_imp_df = pd.DataFrame({'Exch_imp': exch_imp_year['Exch_imp']})
                    exch_exp_df = pd.DataFrame({'Exch_exp': exch_exp_year['Exch_exp']})
                    
                    exch_year = exch_exp_df.merge(exch_imp_df, left_index=True, right_index=True, how='outer')
                    exch_year['Balance'] = exch_year['Exch_exp'] - exch_year['Exch_imp']
                    exch_year = exch_year.assign(Years=year).reset_index()
                    exch_all_years.append(exch_year)
                
                exch_all = pd.concat(exch_all_years)
                exch_all = exch_all.set_index(['Years', 'Resources', 'From', 'To', 'Typical_days', 'Hours']).sort_index()
                self.hourly_results['Exchanges'] = exch_all
            except Exception as e:
                logging.warning(f"Could not save hourly Exchanges: {e}")
        
        return
    
    def get_assets_pathway(self, save_hourly:list=[]):
        """Gets assets (installed capacity and operation) with YEARS dimension
        
        DEFINITIVO: Fixed Storage_tech index finding with dynamic search
        """
        logging.info('Getting Assets (pathway)')
        
        # Get installed capacity F
        f = self.esom.get_var('F')
        
        # Rename indices
        if 'YEARS' in f.index.names:
            f.index = f.index.set_names('Years', level='YEARS')
        if 'REGIONS' in f.index.names:
            f.index = f.index.set_names('Regions', level='REGIONS')
        if 'I in TECHNOLOGIES' in f.index.names:
            f.index = f.index.set_names('Technologies', level='I in TECHNOLOGIES')
        
        # Get F_t (hourly operation)
        f_t = self.esom.get_var('F_t')
        
        # Rename indices
        if 'YEARS' in f_t.index.names:
            f_t.index = f_t.index.set_names('Years', level='YEARS')
        if 'REGIONS' in f_t.index.names:
            f_t.index = f_t.index.set_names('Regions', level='REGIONS')
        if 'I in TECHNOLOGIES' in f_t.index.names:
            f_t.index = f_t.index.set_names('Technologies', level='I in TECHNOLOGIES')
        if 'HOURS' in f_t.index.names:
            f_t.index = f_t.index.set_names('Hours', level='HOURS')
        if 'TYPICAL_DAYS' in f_t.index.names:
            f_t.index = f_t.index.set_names('Typical_days', level='TYPICAL_DAYS')
        
        # Get years and regions
        years_list = sorted(f.index.get_level_values('Years').unique())
        if hasattr(self, 'regions_names'):
            regions_names = self.regions_names
        else:
            regions_names = sorted(f.index.get_level_values('Regions').unique())
        
        # Calculate F_year using temporal expansion (like original esmc.py)
        f_year_list = []

        for year in years_list:
            # Extract F_t for this year
            f_t_year = f_t.xs(year, level='Years')
            
            # Use from_td_to_year to expand to 8760 hours, then sum
            if self.ta is not None:
                f_year_data = self.ta.from_td_to_year(
                    ts_td=f_t_year.reset_index().set_index(['Typical_days', 'Hours'])
                ).groupby(['Regions', 'Technologies']).sum().rename(columns={'F_t': 'F_year'})
            else:
                logging.warning("No temporal aggregation available, using direct TD sum")
                f_year_data = f_t_year.groupby(['Regions', 'Technologies']).sum().rename(columns={'F_t': 'F_year'})
            
            # Add year dimension
            f_year_data = f_year_data.reset_index()
            f_year_data['Years'] = year
            f_year_data = f_year_data.set_index(['Years', 'Regions', 'Technologies'])
            f_year_list.append(f_year_data)

        # Combine all years
        f_year_all = pd.concat(f_year_list)
        
        # Merge F (installed capacity) with F_year (annual operation)
        assets = f.merge(f_year_all, left_index=True, right_index=True, how='left')
        
        # ===== STORAGE ASSETS CALCULATIONS =====
        try:
            # Get Storage_in and Storage_out
            storage_in = self.esom.get_var('Storage_in')
            storage_out = self.esom.get_var('Storage_out')
            
            # Log the actual index names we received from AMPL
            logging.info(f"Storage_in index names from AMPL: {storage_in.index.names}")
            logging.info(f"Storage_out index names from AMPL: {storage_out.index.names}")
            
            # Build a mapping of possible index names to standardized names
            def standardize_index_names(df):
                """Standardize index names to a consistent format"""
                name_mapping = {}
                for idx_name in df.index.names:
                    idx_upper = str(idx_name).upper()
                    # Years
                    if idx_upper == 'YEARS':
                        name_mapping[idx_name] = 'Years'
                    # Regions  
                    elif idx_upper == 'REGIONS':
                        name_mapping[idx_name] = 'Regions'
                    # Storage tech - handle multiple formats
                    elif 'STORAGE' in idx_upper and 'TECH' in idx_upper:
                        name_mapping[idx_name] = 'Storage_tech'
                    # Hours
                    elif idx_upper == 'HOURS':
                        name_mapping[idx_name] = 'Hours'
                    # Typical days - handle multiple formats
                    elif 'TYPICAL' in idx_upper and 'DAY' in idx_upper:
                        name_mapping[idx_name] = 'Typical_days'
                    # Layers
                    elif idx_upper == 'LAYERS':
                        name_mapping[idx_name] = 'Layers'
                
                # Apply the mapping
                if name_mapping:
                    return df.rename_axis(index=name_mapping)
                return df
            
            # Standardize index names
            storage_in = standardize_index_names(storage_in)
            storage_out = standardize_index_names(storage_out)
            
            logging.info(f"After standardization - Storage_in indices: {storage_in.index.names}")
            
            # Check what indices we have
            # For storage: we want to sum over Layers, so don't include it in groupby
            required_indices = ['Years', 'Regions', 'Storage_tech', 'Hours']
            optional_indices = ['Typical_days']  # Typical_days only - NOT Layers (we sum over Layers)
            
            # Build the groupby list dynamically
            groupby_cols = []
            for idx in required_indices:
                if idx in storage_in.index.names:
                    groupby_cols.append(idx)
                else:
                    raise ValueError(f"Required index '{idx}' not found in storage data. Available: {storage_in.index.names}")
            
            # Add optional indices if they exist
            for idx in optional_indices:
                if idx in storage_in.index.names:
                    groupby_cols.append(idx)
            
            logging.info(f"Grouping storage data by: {groupby_cols} (summing over Layers)")
            
            # Group by the identified columns - this will sum over any indices not in groupby_cols (like Layers)
            storage_in_grouped = storage_in.groupby(groupby_cols).sum()
            storage_out_grouped = storage_out.groupby(groupby_cols).sum()
            
            # Process storage for each year
            sto_assets_list = []
            
            for year in years_list:
                # Extract storage data for this year
                storage_in_year = storage_in_grouped.xs(year, level='Years')
                storage_out_year = storage_out_grouped.xs(year, level='Years')
                
                # Compute the balance
                storage_in_df = pd.DataFrame({'Storage_in': storage_in_year['Storage_in']})
                storage_out_df = pd.DataFrame({'Storage_out': storage_out_year['Storage_out']})
                storage_power = storage_out_df.merge(storage_in_df, left_index=True, right_index=True, how='outer')
                storage_power['Storage_power'] = storage_power['Storage_out'].fillna(0) - storage_power['Storage_in'].fillna(0)
                
                # Losses are the sum of the balance over the year
                if self.ta is not None:
                    sto_losses = self.ta.from_td_to_year(
                        ts_td=storage_power[['Storage_power']].reset_index().set_index(['Typical_days', 'Hours'])
                    ).groupby(['Regions', 'Storage_tech']).sum().rename(columns={'Storage_power': 'Losses'})
                else:
                    logging.warning("No temporal aggregation available, using direct TD sum")
                    sto_losses = storage_power[['Storage_power']].groupby(['Regions', 'Storage_tech']).sum().rename(columns={'Storage_power': 'Losses'})

                # Replace Storage_in and Storage_out by values deduced from Storage_power
                threshold = 1e-2
                storage_power['Storage_in'] = storage_power['Storage_power'].mask((storage_power['Storage_power'] > -threshold), np.nan)
                storage_power['Storage_out'] = storage_power['Storage_power'].mask((storage_power['Storage_power'] < threshold), np.nan)

                # Compute total over the year
                if self.ta is not None:
                    sto_flux_year = self.ta.from_td_to_year(
                        ts_td=storage_power[['Storage_out']].reset_index().set_index(['Typical_days', 'Hours'])
                    ).groupby(['Regions', 'Storage_tech']).sum().rename(columns={'Storage_out': 'Year_energy_flux'})
                else:
                    logging.warning("No temporal aggregation available, using direct TD sum")
                    sto_flux_year = storage_power[['Storage_out']].groupby(['Regions', 'Storage_tech']).sum().rename(columns={'Storage_out': 'Year_energy_flux'})
                sto_flux_year.columns = ['Year_energy_flux']
                
                # Merge losses and flux
                sto_year_data = sto_losses.merge(sto_flux_year, left_index=True, right_index=True, how='outer')
                
                # Get F values for storage techs from assets
                f_year = f.xs(year, level='Years')
                f_year_reset = f_year.reset_index()
                f_year_reset = f_year_reset.rename(columns={'Technologies': 'Storage_tech'})
                
                # Merge with storage data
                sto_year_data_reset = sto_year_data.reset_index()
                sto_year_data_reset = sto_year_data_reset.merge(
                    f_year_reset[['Regions', 'Storage_tech', 'F']], 
                    on=['Regions', 'Storage_tech'], 
                    how='left'
                )
                
                # Add year dimension
                sto_year_data_reset['Years'] = year
                sto_year_data_reset = sto_year_data_reset.set_index(['Years', 'Regions', 'Storage_tech'])
                
                sto_assets_list.append(sto_year_data_reset)
            
            # Combine all years
            if sto_assets_list:
                sto_assets = pd.concat(sto_assets_list)
                
                # Rename Storage_tech to Technologies for consistency
                sto_assets = sto_assets.reset_index()
                sto_assets = sto_assets.rename(columns={'Storage_tech': 'Technologies'})
                
                # Set categorical data for sorting
                sto_assets['Years'] = pd.Categorical(sto_assets['Years'], years_list)
                sto_assets['Regions'] = pd.Categorical(sto_assets['Regions'], regions_names)
                tech_order = sorted(sto_assets['Technologies'].unique())
                sto_assets['Technologies'] = pd.Categorical(sto_assets['Technologies'], tech_order)
                
                sto_assets.sort_values(by=['Years', 'Regions', 'Technologies'], axis=0, ignore_index=True, inplace=True)
                sto_assets.set_index(['Years', 'Regions', 'Technologies'], inplace=True)
                
                # Put very small values as nan
                threshold = 1
                # Select only numeric columns for masking
                numeric_cols = sto_assets.select_dtypes(include=[np.number]).columns
                if len(numeric_cols) > 0:
                    sto_assets[numeric_cols] = sto_assets[numeric_cols].mask(
                        (sto_assets[numeric_cols] > -threshold) & (sto_assets[numeric_cols] < threshold), 
                        np.nan
                    )
                
                # Update F_year in assets df for STORAGE_TECH
                for idx in sto_assets.index:
                    if idx in assets.index:
                        assets.loc[idx, 'F_year'] = sto_assets.loc[idx, 'Losses']
                
                logging.info(f"Successfully calculated storage assets")
            else:
                sto_assets = pd.DataFrame()
                
        except Exception as e:
            logging.error(f"Could not calculate storage assets: {e}")
            import traceback
            traceback.print_exc()
            sto_assets = pd.DataFrame()
        
        # Set categorical data for sorting assets
        assets = assets.reset_index()
        assets['Years'] = pd.Categorical(assets['Years'], years_list)
        assets['Regions'] = pd.Categorical(assets['Regions'], regions_names)
        
        tech_order = sorted(assets['Technologies'].unique())
        assets['Technologies'] = pd.Categorical(assets['Technologies'], tech_order)
        
        assets.sort_values(by=['Years', 'Regions', 'Technologies'], axis=0, ignore_index=True, inplace=True)
        assets.set_index(['Years', 'Regions', 'Technologies'], inplace=True)
        
        # Put very small values as nan
        threshold = 1e-5
        # Select only numeric columns for masking
        numeric_cols = assets.select_dtypes(include=[np.number]).columns
        if len(numeric_cols) > 0:
            assets[numeric_cols] = assets[numeric_cols].mask(
                (assets[numeric_cols] > -threshold) & (assets[numeric_cols] < threshold), 
                np.nan
            )
        
        # Store results
        self.results['Assets'] = assets
        self.results_all['Assets'] = assets.groupby(['Years', 'Technologies'], observed=True).sum()
        
        self.results['Sto_assets'] = sto_assets
        if not sto_assets.empty:
            self.results_all['Sto_assets'] = sto_assets.groupby(['Years', 'Technologies'], observed=True).sum()
        else:
            self.results_all['Sto_assets'] = pd.DataFrame()
        
        # Save hourly if requested
        if 'Assets' in save_hourly:
            try:
                self.hourly_results['Assets'] = f_t.reset_index().set_index(
                    ['Years', 'Regions', 'Technologies', 'Typical_days', 'Hours']).sort_index()
            except Exception as e:
                logging.warning(f"Could not save hourly Assets: {e}")
        
        if 'Storage' in save_hourly:
            try:
                # Only try to save if storage variables were successfully processed
                if 'storage_in_grouped' in locals() and 'storage_out_grouped' in locals():
                    # Save storage power hourly data
                    storage_power_all_years = []
                    for year in years_list:
                        storage_in_year = storage_in_grouped.xs(year, level='Years')
                        storage_out_year = storage_out_grouped.xs(year, level='Years')
                        
                        storage_in_df = pd.DataFrame({'Storage_in': storage_in_year['Storage_in']})
                        storage_out_df = pd.DataFrame({'Storage_out': storage_out_year['Storage_out']})
                        storage_power_year = storage_out_df.merge(storage_in_df, left_index=True, right_index=True, how='outer')
                        storage_power_year['Storage_power'] = storage_power_year['Storage_out'] - storage_power_year['Storage_in']
                        storage_power_year = storage_power_year.assign(Years=year).reset_index()
                        storage_power_all_years.append(storage_power_year)
                    
                    storage_power_all = pd.concat(storage_power_all_years)
                    storage_power_all = storage_power_all.set_index(['Years', 'Regions', 'Storage_tech', 'Typical_days', 'Hours']).sort_index()
                    self.hourly_results['Storage_power'] = storage_power_all
                    
                    # Also try to save Storage_level
                    try:
                        storage_level = self.esom.get_var('Storage_level')
                        if 'YEARS' in storage_level.index.names:
                            storage_level.index = storage_level.index.set_names('Years', level='YEARS')
                        self.hourly_results['Storage_level'] = storage_level
                    except:
                        pass
            except Exception as e:
                logging.warning(f"Could not save hourly Storage: {e}")
        
        return
    
    def get_year_balance_pathway(self):
        """Gets year balance (energy balance) with YEARS dimension
        
        DEFINITIVO: Completely reimplemented to match the one-year model structure:
        - Rows: Technologies + Resources + END_USES
        - Columns: Layers (ELECTRICITY, GAS, etc.)
        - Values: Energy balance for each layer
        
        Matches esmc.py get_year_balance() but with Years dimension added
        """
        logging.info('Getting Year_balance (pathway)')
        
        # Get years and regions lists
        if hasattr(self, 'years'):
            years_list = self.years
        else:
            try:
                years_list = sorted(self.results['Assets'].index.get_level_values('Years').unique())
            except:
                years_list = []
        
        if hasattr(self, 'regions_names'):
            regions_names = self.regions_names
        else:
            try:
                regions_names = sorted(self.results['Assets'].index.get_level_values('Regions').unique())
            except:
                regions_names = []
        
        # Initialize list to collect year_balance for each year
        year_balance_list = []
        
        for year in years_list:
            logging.info(f"Processing year_balance for year: {year}")
            
            # ===== EXTRACT END_USES =====
            try:
                end_uses_td = self.esom.get_var('End_uses')
                
                # Rename indices
                if 'YEARS' in end_uses_td.index.names:
                    end_uses_td.index = end_uses_td.index.set_names('Years', level='YEARS')
                if 'REGIONS' in end_uses_td.index.names:
                    end_uses_td.index = end_uses_td.index.set_names('Regions', level='REGIONS')
                if 'LAYERS' in end_uses_td.index.names:
                    end_uses_td.index = end_uses_td.index.set_names('Layers', level='LAYERS')
                if 'HOURS' in end_uses_td.index.names:
                    end_uses_td.index = end_uses_td.index.set_names('Hours', level='HOURS')
                if 'TYPICAL_DAYS' in end_uses_td.index.names:
                    end_uses_td.index = end_uses_td.index.set_names('Typical_days', level='TYPICAL_DAYS')
                
                # Extract for this year - handle type mismatch and YEAR_ prefix
                years_in_index = end_uses_td.index.get_level_values('Years').unique()
                year_to_use = year
                
                # Match the year format used in the index
                if len(years_in_index) > 0:
                    sample_year = years_in_index[0]
                    
                    # Handle YEAR_ prefix format (e.g., 'YEAR_2015')
                    if isinstance(sample_year, str) and sample_year.startswith('YEAR_'):
                        # Index has YEAR_ prefix
                        if isinstance(year, str) and year.startswith('YEAR_'):
                            year_to_use = year  # Already has prefix
                        else:
                            # Add YEAR_ prefix
                            year_num = str(year).replace('YEAR_', '')  # Remove any existing prefix
                            year_to_use = f'YEAR_{year_num}'
                    # Handle plain string vs int
                    elif isinstance(sample_year, str) and not isinstance(year, str):
                        year_to_use = str(year)
                    elif isinstance(sample_year, (int, np.integer)) and isinstance(year, str):
                        year_to_use = int(year.replace('YEAR_', ''))  # Remove prefix if present
                
                end_uses_year = end_uses_td.xs(year_to_use, level='Years')
                
                # Convert to yearly using temporal expansion
                if self.ta is not None:
                    end_uses_year_data = self.ta.from_td_to_year(
                        ts_td=end_uses_year.reset_index().set_index(['Typical_days', 'Hours'])
                    ).groupby(['Regions', 'Layers']).sum().rename(columns={'End_uses': 'End_uses'})
                else:
                    logging.warning("No temporal aggregation available, using direct TD sum")
                    end_uses_year_data = end_uses_year.groupby(['Regions', 'Layers']).sum()

                # Negate for balance (end uses are consumption)
                end_uses_year_data = -end_uses_year_data
                
                # Pivot to have Layers as columns
                end_uses_pivot = end_uses_year_data.reset_index()
                end_uses_pivot['Elements'] = 'END_USES'
                
                # Pivot with Layers as columns
                end_uses_pivot = end_uses_pivot.pivot_table(
                    index=['Regions', 'Elements'],
                    columns='Layers',
                    values='End_uses',
                    aggfunc='sum'
                )
                
            except Exception as e:
                logging.warning(f"Could not get End_uses for year {year}: {e}")
                end_uses_pivot = pd.DataFrame()
            
            # ===== EXTRACT ASSETS AND RESOURCES YEARLY FLUXES =====
            # Get F_year from assets
            if self.results['Assets'] is not None and not self.results['Assets'].empty:
                try:
                    # Match year format to Assets index
                    years_in_assets = self.results['Assets'].index.get_level_values('Years').unique()
                    year_to_use = year
                    if len(years_in_assets) > 0:
                        sample_year = years_in_assets[0]
                        if isinstance(sample_year, str) and sample_year.startswith('YEAR_'):
                            year_num = str(year).replace('YEAR_', '')
                            year_to_use = f'YEAR_{year_num}'
                        elif isinstance(sample_year, str) and not isinstance(year, str):
                            year_to_use = str(year)
                        elif isinstance(sample_year, (int, np.integer)) and isinstance(year, str):
                            year_to_use = int(year.replace('YEAR_', ''))
                    
                    assets_year = self.results['Assets'].xs(year_to_use, level='Years')
                    f_year = assets_year['F_year'].reset_index() \
                        .rename(columns={'Technologies': 'Elements'})
                except Exception as e:
                    logging.warning(f"Could not extract Assets for year {year}: {e}")
                    f_year = pd.DataFrame()
            else:
                f_year = pd.DataFrame()
            
            # Get R_year from resources
            if self.results['Resources'] is not None and not self.results['Resources'].empty:
                try:
                    # Match year format to Resources index
                    years_in_resources = self.results['Resources'].index.get_level_values('Years').unique()
                    year_to_use = year
                    if len(years_in_resources) > 0:
                        sample_year = years_in_resources[0]
                        if isinstance(sample_year, str) and sample_year.startswith('YEAR_'):
                            year_num = str(year).replace('YEAR_', '')
                            year_to_use = f'YEAR_{year_num}'
                        elif isinstance(sample_year, str) and not isinstance(year, str):
                            year_to_use = str(year)
                        elif isinstance(sample_year, (int, np.integer)) and isinstance(year, str):
                            year_to_use = int(year.replace('YEAR_', ''))
                    
                    resources_year = self.results['Resources'].xs(year_to_use, level='Years')
                    
                    # Calculate total resource usage (local + exterior + import - export)
                    r_year_total = pd.Series(0.0, index=resources_year.index)
                    if 'R_year_local' in resources_year.columns:
                        r_year_total += resources_year['R_year_local'].fillna(0)
                    if 'R_year_exterior' in resources_year.columns:
                        r_year_total += resources_year['R_year_exterior'].fillna(0)
                    if 'R_year_import' in resources_year.columns:
                        r_year_total += resources_year['R_year_import'].fillna(0)
                    if 'R_year_export' in resources_year.columns:
                        r_year_total -= resources_year['R_year_export'].fillna(0)
                    
                    r_year = r_year_total.reset_index().rename(columns={'Resources': 'Elements', 0: 'R_year'})
                except:
                    r_year = pd.DataFrame()
            else:
                r_year = pd.DataFrame()
            
            # Combine F_year and R_year
            if not f_year.empty and not r_year.empty:
                year_fluxes = pd.concat([
                    r_year.set_index(['Regions', 'Elements'])['R_year'],
                    f_year.set_index(['Regions', 'Elements'])['F_year']
                ], axis=0)
            elif not f_year.empty:
                year_fluxes = f_year.set_index(['Regions', 'Elements'])['F_year']
            elif not r_year.empty:
                year_fluxes = r_year.set_index(['Regions', 'Elements'])['R_year']
            else:
                year_fluxes = pd.Series(dtype=float)
            
            # ===== GET LAYERS_IN_OUT PARAMETER =====
            try:
                # Get layers_in_out parameter from AMPL for this year
                layers_in_out_param = self.esom.get_param('layers_in_out')
                
                # Filter for this year if it has Years dimension
                if 'YEARS' in layers_in_out_param.index.names or any('YEAR' in str(name).upper() for name in layers_in_out_param.index.names):
                    # Rename Years index
                    for name in layers_in_out_param.index.names:
                        if 'YEAR' in str(name).upper():
                            layers_in_out_param.index = layers_in_out_param.index.set_names('Years', level=name)
                            break
                    
                    # Extract for this year - handle type mismatch and YEAR_ prefix
                    years_in_index = layers_in_out_param.index.get_level_values('Years').unique()
                    year_to_use = year
                    
                    # Match the year format used in the index
                    if len(years_in_index) > 0:
                        sample_year = years_in_index[0]
                        
                        # Handle YEAR_ prefix format (e.g., 'YEAR_2015')
                        if isinstance(sample_year, str) and sample_year.startswith('YEAR_'):
                            # Index has YEAR_ prefix
                            if isinstance(year, str) and year.startswith('YEAR_'):
                                year_to_use = year  # Already has prefix
                            else:
                                # Add YEAR_ prefix
                                year_num = str(year).replace('YEAR_', '')  # Remove any existing prefix
                                year_to_use = f'YEAR_{year_num}'
                        # Handle plain string vs int
                        elif isinstance(sample_year, str) and not isinstance(year, str):
                            year_to_use = str(year)
                        elif isinstance(sample_year, (int, np.integer)) and isinstance(year, str):
                            year_to_use = int(year.replace('YEAR_', ''))  # Remove prefix if present
                    
                    layers_in_out_year = layers_in_out_param.xs(year_to_use, level='Years')
                else:
                    # No year dimension, use as-is
                    layers_in_out_year = layers_in_out_param
                
                # Rename indices
                index_names = list(layers_in_out_year.index.names)
                new_names = []
                for name in index_names:
                    if 'TECHNOLOGIES' in str(name).upper() or 'RESOURCES' in str(name).upper():
                        new_names.append('Elements')
                    elif 'LAYERS' in str(name).upper():
                        new_names.append('Layers')
                    else:
                        new_names.append(name)
                layers_in_out_year.index.names = new_names
                
                # Pivot to have Layers as columns
                layers_in_out_pivot = layers_in_out_year.reset_index().pivot_table(
                    index='Elements',
                    columns='Layers',
                    values='layers_in_out',
                    aggfunc='first'
                )
                
            except Exception as e:
                logging.warning(f"Could not get layers_in_out for year {year}: {e}")
                layers_in_out_pivot = pd.DataFrame()
            
            # ===== CALCULATE YEAR BALANCE =====
            if not layers_in_out_pivot.empty and not year_fluxes.empty:
                # For each region, multiply year_fluxes by layers_in_out
                year_balance_parts = []
                
                for region in regions_names:
                    # Get fluxes for this region
                    if region in year_fluxes.index.get_level_values('Regions'):
                        region_fluxes = year_fluxes.xs(region, level='Regions')
                        
                        # Multiply by layers_in_out
                        region_balance = layers_in_out_pivot.mul(region_fluxes, axis=0)
                        
                        # Add region index
                        region_balance['Regions'] = region
                        region_balance = region_balance.reset_index().set_index(['Regions', 'Elements'])
                        
                        year_balance_parts.append(region_balance)
                
                if year_balance_parts:
                    tech_res_balance = pd.concat(year_balance_parts)
                else:
                    tech_res_balance = pd.DataFrame()
            else:
                tech_res_balance = pd.DataFrame()
            
            # ===== COMBINE TECHNOLOGIES/RESOURCES WITH END_USES =====
            if not tech_res_balance.empty and not end_uses_pivot.empty:
                # Make sure columns match
                all_layers = sorted(set(tech_res_balance.columns) | set(end_uses_pivot.columns))
                
                # Reindex both to have same columns
                tech_res_balance = tech_res_balance.reindex(columns=all_layers)
                end_uses_pivot = end_uses_pivot.reindex(columns=all_layers)
                
                # Concatenate
                year_balance_year = pd.concat([tech_res_balance, end_uses_pivot], axis=0)
            elif not tech_res_balance.empty:
                year_balance_year = tech_res_balance
            elif not end_uses_pivot.empty:
                year_balance_year = end_uses_pivot
            else:
                year_balance_year = pd.DataFrame()
            
            # Add year dimension
            if not year_balance_year.empty:
                year_balance_year = year_balance_year.assign(Years=year).reset_index()
                year_balance_year = year_balance_year.set_index(['Years', 'Regions', 'Elements'])
                year_balance_list.append(year_balance_year)
        
        # ===== COMBINE ALL YEARS =====
        if year_balance_list:
            year_balance = pd.concat(year_balance_list)
            
            # Set categorical data for sorting
            year_balance = year_balance.reset_index()
            year_balance['Years'] = pd.Categorical(year_balance['Years'], years_list)
            year_balance['Regions'] = pd.Categorical(year_balance['Regions'], regions_names)
            
            # Get elements order (technologies, resources, END_USES)
            elements_order = sorted([e for e in year_balance['Elements'].unique() if e != 'END_USES'])
            elements_order.append('END_USES')  # END_USES should be last
            year_balance['Elements'] = pd.Categorical(year_balance['Elements'], elements_order)
            
            year_balance.sort_values(by=['Years', 'Regions', 'Elements'], axis=0, ignore_index=True, inplace=True)
            year_balance.set_index(['Years', 'Regions', 'Elements'], inplace=True)
            
            # Put very small values as nan
            threshold = 1e-5
            year_balance = year_balance.mask(
                (year_balance.min(axis=1) > -threshold) & (year_balance.max(axis=1) < threshold),
                np.nan
            )
            
            logging.info(f"Successfully calculated Year_balance with {len(year_balance)} rows")
        else:
            year_balance = pd.DataFrame()
            logging.warning("Year_balance is empty")
        
        # Store results
        self.results['Year_balance'] = year_balance
        if not year_balance.empty:
            self.results_all['Year_balance'] = year_balance.groupby(['Years', 'Elements'], observed=True).sum()
        else:
            self.results_all['Year_balance'] = pd.DataFrame()
        
        return
    
    
    def get_curt_pathway(self, save_hourly:list=[]):
        """Gets yearly curtailment of renewables with YEARS dimension"""
        logging.info('Getting Curt (pathway)')
        
        try:
            curt_t = self.esom.get_var('Curt')
            
            if 'YEARS' in curt_t.index.names:
                curt_t.index = curt_t.index.set_names('Years', level='YEARS')
            if 'REGIONS' in curt_t.index.names:
                curt_t.index = curt_t.index.set_names('Regions', level='REGIONS')
            if 'I in TECHNOLOGIES' in curt_t.index.names:
                curt_t.index = curt_t.index.set_names('Technologies', level='I in TECHNOLOGIES')
            if 'HOURS' in curt_t.index.names:
                curt_t.index = curt_t.index.set_names('Hours', level='HOURS')
            if 'TYPICAL_DAYS' in curt_t.index.names:
                curt_t.index = curt_t.index.set_names('Typical_days', level='TYPICAL_DAYS')
            
            years_list = sorted(curt_t.index.get_level_values('Years').unique())
            if hasattr(self, 'regions_names'):
                regions_names = self.regions_names
            else:
                regions_names = sorted(curt_t.index.get_level_values('Regions').unique())
            
            # Calculate yearly curtailment using t_op weighting
            curt_list = []
            
            for year in years_list:
                curt_t_year = curt_t.xs(year, level='Years')
                
                # Use from_td_to_year to convert to yearly values
                if self.ta is not None:
                    curt_year = self.ta.from_td_to_year(
                        ts_td=curt_t_year.reset_index().set_index(['Typical_days', 'Hours'])
                    ).groupby(['Regions', 'Technologies']).sum().rename(columns={'Curt': 'Curt'})
                else:
                    logging.warning("No temporal aggregation available, using direct TD sum")
                    curt_year = curt_t_year.groupby(['Regions', 'Technologies']).sum()

                # Add year dimension
                curt_year = curt_year.reset_index()
                curt_year['Years'] = year
                curt_year = curt_year.set_index(['Years', 'Regions', 'Technologies'])
                curt_list.append(curt_year)
            
            curt = pd.concat(curt_list)
            
            curt.reset_index(inplace=True)
            curt['Years'] = pd.Categorical(curt['Years'], years_list)
            curt['Regions'] = pd.Categorical(curt['Regions'], regions_names)
            
            tech_order = sorted(curt['Technologies'].unique())
            curt['Technologies'] = pd.Categorical(curt['Technologies'], tech_order)
            
            curt.set_index(['Years', 'Regions', 'Technologies'], inplace=True)
            curt.sort_index(inplace=True)
            
            self.results['Curt'] = curt
            self.results_all['Curt'] = curt.groupby(['Years', 'Technologies'], observed=True).sum()
            
            if 'Curt' in save_hourly:
                try:
                    self.hourly_results['Curt'] = curt_t.reset_index().set_index(
                        ['Years', 'Regions', 'Technologies', 'Typical_days', 'Hours']).sort_index()
                except Exception as e:
                    logging.warning(f"Could not save hourly Curt: {e}")
            
        except Exception as e:
            logging.warning(f"Curtailment variable not available or error: {e}")
            self.results['Curt'] = pd.DataFrame()
            self.results_all['Curt'] = pd.DataFrame()
        
        return    
    
    def collect_new_old_decom(self, output_dir):
        """
        Collect F_new, F_old, F_decom from AMPL and save to CSV
        
        Parameters
        ----------
        output_dir : Path
            Directory to save results
            
        Returns
        -------
        pd.DataFrame
            Combined dataframe with F_new, F_old, F_decom
        """
        print("\n" + "="*70)
        print("COLLECTING NEW/OLD/DECOM DATA")
        print("="*70)
        
        try:
            ampl = self.esom.ampl
            regions_names = self.regions_names
            
            # Get F_new
            print("Extracting F_new...")
            f_new_ampl = ampl.getVariable('F_new')
            f_new_df = f_new_ampl.getValues().toPandas()
            
            print(f"  F_new raw shape: {f_new_df.shape}")
            print(f"  F_new index names: {f_new_df.index.names}")
            print(f"  F_new columns: {list(f_new_df.columns)}")
            
            # The DataFrame already has MultiIndex set
            # Just rename the column and index names
            value_col = [c for c in f_new_df.columns if '.val' in str(c)][0]
            f_new_df = f_new_df[[value_col]].copy()
            f_new_df.columns = ['F_new']
            
            # Determine if multi-regional based on number of index levels
            n_levels = f_new_df.index.nlevels
            
            if n_levels == 3:
                # Multi-regional: (Phase, Region, Technology)
                f_new_df.index.names = ['Phases', 'Regions', 'Technologies']
                is_multiregional = True
                print(f"  Model type: Multi-regional (3 index levels)")
            elif n_levels == 2:
                # Single-regional: (Phase, Technology)
                f_new_df.index.names = ['Phases', 'Technologies']
                is_multiregional = False
                print(f"  Model type: Single-regional (2 index levels)")
            else:
                print(f"  WARNING: Unexpected number of index levels: {n_levels}")
                is_multiregional = n_levels > 2
            
            print(f"  F_new processed: {f_new_df.shape} - {f_new_df.index.names}")
            
            # Get F_old
            print("\nExtracting F_old...")
            f_old_ampl = ampl.getVariable('F_old')
            f_old_df = f_old_ampl.getValues().toPandas()
            
            value_col = [c for c in f_old_df.columns if '.val' in str(c)][0]
            f_old_df = f_old_df[[value_col]].copy()
            f_old_df.columns = ['F_old']
            f_old_df.index.names = f_new_df.index.names  # Use same names as F_new
            
            print(f"  F_old processed: {f_old_df.shape} - {f_old_df.index.names}")
            
            # Get F_decom
            print("\nExtracting F_decom...")
            f_decom_ampl = ampl.getVariable('F_decom')
            f_decom_df = f_decom_ampl.getValues().toPandas()
            
            print(f"  F_decom raw shape: {f_decom_df.shape}")
            print(f"  F_decom index levels: {f_decom_df.index.nlevels}")
            print(f"  F_decom index names: {f_decom_df.index.names}")
            
            value_col = [c for c in f_decom_df.columns if '.val' in str(c)][0]
            f_decom_df = f_decom_df[[value_col]].copy()
            f_decom_df.columns = ['F_decom']
            
            print(f"  F_decom before aggregation: {f_decom_df.shape}, levels: {f_decom_df.index.nlevels}")
            
            # Aggregate F_decom over phase_built dimension
            # F_decom has structure: (phase_decom, phase_built, region, tech) or (phase_decom, phase_built, tech)
            if is_multiregional:
                if f_decom_df.index.nlevels == 4:
                    # (phase_decom, phase_built, region, tech) -> (phase_decom, region, tech)
                    print("  Aggregating F_decom: 4 levels -> 3 levels")
                    f_decom_df = f_decom_df.groupby(level=[0, 2, 3]).sum()
                    f_decom_df.index.names = ['Phases', 'Regions', 'Technologies']
                elif f_decom_df.index.nlevels == 3:
                    print("  F_decom already has 3 levels, assuming correct structure")
                    f_decom_df.index.names = ['Phases', 'Regions', 'Technologies']
                else:
                    print(f"  WARNING: Unexpected F_decom levels: {f_decom_df.index.nlevels}")
                    f_decom_df.index.names = f_new_df.index.names
            else:
                if f_decom_df.index.nlevels == 3:
                    # (phase_decom, phase_built, tech) -> (phase_decom, tech)
                    print("  Aggregating F_decom: 3 levels -> 2 levels")
                    f_decom_df = f_decom_df.groupby(level=[0, 2]).sum()
                    f_decom_df.index.names = ['Phases', 'Technologies']
                elif f_decom_df.index.nlevels == 2:
                    print("  F_decom already has 2 levels, assuming correct structure")
                    f_decom_df.index.names = ['Phases', 'Technologies']
                else:
                    print(f"  WARNING: Unexpected F_decom levels: {f_decom_df.index.nlevels}")
                    f_decom_df.index.names = f_new_df.index.names
            
            print(f"  F_decom aggregated: {f_decom_df.shape} - {f_decom_df.index.names}")
            
            # Merge all
            print("\nMerging dataframes...")
            new_old_decom = f_new_df.join(f_old_df, how='outer').join(f_decom_df, how='outer')
            new_old_decom.fillna(0, inplace=True)
            
            print(f"  Merged shape: {new_old_decom.shape}")
            print(f"  Merged index: {new_old_decom.index.names}")
            
            # Remove very small values
            threshold = 1e-5
            new_old_decom = new_old_decom.mask((new_old_decom.abs() < threshold), np.nan)
            
            # Get storage technologies and filter them out
            try:
                sto_tech = list(ampl.getSet('STORAGE_TECH').getValues())
                print(f"  Filtering {len(sto_tech)} storage technologies")
                new_old_decom = new_old_decom.loc[
                    ~new_old_decom.index.get_level_values('Technologies').isin(sto_tech), :
                ]
                print(f"  After filtering: {new_old_decom.shape}")
            except Exception as e:
                print(f"  Warning: Could not filter storage technologies: {e}")
            
            print(f"\nFinal shape: {new_old_decom.shape}")
            print(f"Index: {new_old_decom.index.names}")
            print(f"Non-null values:")
            print(f"  F_new: {new_old_decom['F_new'].notna().sum()}")
            print(f"  F_old: {new_old_decom['F_old'].notna().sum()}")
            print(f"  F_decom: {new_old_decom['F_decom'].notna().sum()}")
            
            # Show sample data
            print(f"\nSample data (first 10 non-zero rows):")
            sample = new_old_decom[(new_old_decom['F_new'].notna()) | 
                                   (new_old_decom['F_old'].notna()) | 
                                   (new_old_decom['F_decom'].notna())].head(10)
            print(sample)
            
            # Save to CSV
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            
            csv_path = output_dir / 'New_old_decom.csv'
            new_old_decom.to_csv(csv_path, sep=';')
            print(f"\nâœ“ Saved to: {csv_path}")
            
            return new_old_decom
            
        except Exception as e:
            print(f"âœ— Error collecting data: {e}")
            import traceback
            traceback.print_exc()
            return None
    
    
    @staticmethod
    def _get_sector_mappings():
        """
        Get technology to sector mappings
        
        Returns
        -------
        tuple
            (tech_to_end_use dict, end_use_to_sector dict)
        """
        tech_by_end_use = {
            "ELECTRICITY": ["NUCLEAR", "CCGT", "CCGT_AMMONIA", "COAL_US", "COAL_IGCC", 
                            "PV_UTILITY", "WIND_ONSHORE", "WIND_OFFSHORE", "HYDRO_DAM", "HYDRO_RIVER", 
                            "GEOTHERMAL", "OCGT", "GENSET_DIESEL", "ST_BIOMASS", "ST_SNG",
                            "FUEL_CELL", "FB_ST_BIOMASS", "CFB_ST_BIOMASS", "BFB_ST_BIOMASS"],
            "LIGHTING_R_C": ["CONVENTIONAL_BULB", "LED_BULB"],
            "LIGHTING_P": ["CONVENTIONAL_LIGHT", "LED_LIGHT"],
            "HEAT_HIGH_T": ["IND_COGEN_GAS", "IND_COGEN_WOOD", "IND_COGEN_WASTE", 
                            "IND_BOILER_GAS", "IND_BOILER_WOOD", "IND_BOILER_OIL", 
                            "IND_BOILER_COAL", "IND_BOILER_WASTE", "IND_DIRECT_ELEC",
                            "IND_BOILER_DIESEL", "IND_BOILER_LPG"],
            "HEAT_LOW_T_DHN": ["DHN_HP_ELEC", "DHN_COGEN_GAS", "DHN_COGEN_WOOD", 
                               "DHN_COGEN_WET_BIOMASS", "DHN_COGEN_BIO_HYDROLYSIS", 
                               "DHN_COGEN_WASTE", "DHN_BOILER_GAS", "DHN_BOILER_WOOD",
                               "DHN_BOILER_OIL", "DHN_DEEP_GEO", "DHN_SOLAR"],
            "HEAT_LOW_T_DECEN": ["DEC_HP_ELEC", "DEC_THHP_GAS", "DEC_COGEN_GAS", 
                                 "DEC_COGEN_OIL", "DEC_ADVCOGEN_GAS", "DEC_ADVCOGEN_H2", 
                                 "DEC_BOILER_GAS", "DEC_BOILER_WOOD", "DEC_BOILER_OIL",
                                 "DEC_SOLAR", "DEC_DIRECT_ELEC"],
            "PROCESS_COOLING": ["IND_ELEC_COLD"],
            "SPACE_COOLING": ["DEC_THHP_GAS_COLD", "DEC_ELEC_COLD"],
            "FOOD_PRESERVATION": ["REFRIGERATOR_EL"],
            "COOKING": ["STOVE_WOOD", "STOVE_LPG", "STOVE_NG", "STOVE_OIL", "STOVE_ELEC"],
            "MOB_PUBLIC": ["TRAMWAY_TROLLEY", "BUS_COACH_DIESEL", "BUS_COACH_GASOLINE", 
                           "BUS_COACH_HYDIESEL", "BUS_COACH_CNG_STOICH", 
                           "BUS_COACH_FC_HYBRIDH2", "TRAIN_PUB"],
            "MOB_PRIVATE": ["CAR_GASOLINE", "CAR_DIESEL", "CAR_NG", "CAR_METHANOL", 
                            "CAR_HEV", "CAR_PHEV", "CAR_BEV", "CAR_FUEL_CELL"],
            "MOB_FREIGHT_RAIL": ["TRAIN_FREIGHT", "TRAIN_FREIGHT_ELEC"],
            "MOB_FREIGHT_BOAT": ["BOAT_FREIGHT_DIESEL", "BOAT_FREIGHT_NG", 
                                 "BOAT_FREIGHT_METHANOL", "BOAT_FREIGHT_ELEC"],
            "MOB_FREIGHT_ROAD": ["TRUCK_DIESEL", "TRUCK_FUEL_CELL", "TRUCK_NG", 
                                 "TRUCK_METHANOL", "TRUCK_ELEC", "TRUCK_FG", "VAN_FG"],
            "MOB_FREIGHT_AIR": ["PLANE"],
            "MECHANICAL_ENERGY": ["COMM_MACHINERY_DIESEL", "COMM_MACHINERY_EL",
                                  "IND_MACHINERY_EL", "TRACTOR_DIESEL", "TRACTOR_EL",
                                  "AGR_MACHINERY_DIESEL", "AGR_MACHINERY_EL",
                                  "MIN_MACHINERY_DIESEL", "MIN_MACHINERY_EL",
                                  "FISH_MACHINERY_DIESEL", "FISH_MACHINERY_EL"],
            "NON_ENERGY": ["HABER_BOSCH", "SYN_METHANOLATION", "METHANE_TO_METHANOL",
                           "BIOMASS_TO_METHANOL", "OIL_TO_HVC", "GAS_TO_HVC", 
                           "BIOMASS_TO_HVC", "METHANOL_TO_HVC"],
            "BIOFUELS": ["H2_ELECTROLYSIS", "H2_NG", "H2_BIOMASS", "BIOMASS_TO_METHANE"],
        }
        
        sector_mapping = {
            "ELECTRICITY": ["ELECTRICITY"],
            "LIGHTING": ["LIGHTING_R_C", "LIGHTING_P"],
            "HEAT_HIGH_T": ["HEAT_HIGH_T"],
            "HEAT_LOW_T": ["HEAT_LOW_T_DHN", "HEAT_LOW_T_DECEN"],
            "COOLING": ["PROCESS_COOLING", "SPACE_COOLING"],
            "FOOD_PRESERVATION": ["FOOD_PRESERVATION"],
            "COOKING": ["COOKING"],
            "MOBILITY_PASSENGER": ["MOB_PUBLIC", "MOB_PRIVATE"],
            "MOBILITY_FREIGHT": ["MOB_FREIGHT_RAIL", "MOB_FREIGHT_BOAT", 
                                 "MOB_FREIGHT_ROAD", "MOB_FREIGHT_AIR"],
            "MECHANICAL_ENERGY": ["MECHANICAL_ENERGY"],
            "NON_ENERGY": ["NON_ENERGY"],
            "BIOFUELS": ["BIOFUELS"],
        }
        
        # Create reverse mappings
        tech_to_end_use = {}
        for end_use, techs in tech_by_end_use.items():
            for tech in techs:
                tech_to_end_use[tech] = end_use
        
        end_use_to_sector = {}
        for sector, end_uses in sector_mapping.items():
            for end_use in end_uses:
                end_use_to_sector[end_use] = sector
        
        return tech_to_end_use, end_use_to_sector
    
    
    def plot_new_old_decom_by_sector(self, new_old_decom_df, output_dir):
        """
        Generate plots by sector for new/old/decom data
        
        Parameters
        ----------
        new_old_decom_df : pd.DataFrame
            DataFrame with F_new, F_old, F_decom (indexed by Phases, [Regions], Technologies)
        output_dir : Path
            Directory to save plots
            
        Returns
        -------
        pd.DataFrame
            Full plotting data
        """
        print("\n" + "="*70)
        print("GENERATING NEW/OLD/DECOM PLOTS")
        print("="*70)
        
        try:
            import plotly.express as px
            import plotly.io as pio
        except ImportError:
            print("✗ Plotly not installed. Install with: pip install plotly kaleido")
            print("  Skipping plot generation, but data is saved in CSV")
            return None
        
        pio.renderers.default = 'browser'
        
        # Detect if multi-regional
        is_multiregional = 'Regions' in new_old_decom_df.index.names
        print(f"Model type: {'Multi-regional' if is_multiregional else 'Single-region'}")
        
        # Prepare data for plotting
        F_new = new_old_decom_df[['F_new']].dropna(how='all').reset_index()
        F_new['Type'] = 'F_new'
        F_new.rename(columns={'F_new': 'Tech_Cap'}, inplace=True)
        
        F_old = new_old_decom_df[['F_old']].dropna(how='all').reset_index()
        F_old['Type'] = 'F_old'
        F_old.rename(columns={'F_old': 'Tech_Cap'}, inplace=True)
        
        F_decom = new_old_decom_df[['F_decom']].dropna(how='all').reset_index()
        F_decom['Type'] = 'F_decom'
        F_decom.rename(columns={'F_decom': 'Tech_Cap'}, inplace=True)
        
        # Combine
        df_full = pd.concat([F_new, F_old, F_decom], ignore_index=True)
        
        print(f"Combined data shape: {df_full.shape}")
        print(f"Phases: {sorted(df_full['Phases'].unique())}")
        if is_multiregional:
            print(f"Regions: {sorted(df_full['Regions'].unique())}")
        print(f"Technologies: {df_full['Technologies'].nunique()} unique")
        
        # Map technologies to categories using the same mapping as cost functions
        tech_to_category = self._get_tech_to_category_mapping()
        df_full['Category'] = df_full['Technologies'].map(tech_to_category)
        
        # Apply display name mapping (same as cost functions)
        category_display_map = {
            'ELECTRICITY': 'ELECTRICITY',
            'LIGHTING': 'LIGHTING',
            'HEAT_HIGH_T': 'HEAT HIGH TEMPERATURE',
            'HEAT_LOW_T': 'HEAT LOW TEMPERATURE',
            'COOLING': 'COOLING & REFRIGERATION',
            'FOOD_PRESERVATION': 'COOLING & REFRIGERATION',
            'COOKING': 'COOKING',
            'MOBILITY_PASSENGER': 'MOBILITY PASSENGER',
            'MOBILITY_FREIGHT': 'MOBILITY FREIGHT',
            'MECHANICAL_ENERGY': 'MECHANICAL ENERGY',
            'NON_ENERGY': 'NON-ENERGY',
            'BIOFUELS_EFUELS': 'BIOFUELS & EFUELS PRODUCTION',
            'ELEC_INFRASTRUCTURE': 'ELECTRICITY INFRASTRUCTURE',
            'HYDROCARBON_INFRASTRUCTURE': 'HYDROCARBON INFRASTRUCTURE',
            'H2_INFRASTRUCTURE': 'H2 INFRASTRUCTURE',
            'OTHER_INFRASTRUCTURE': 'OTHER INFRASTRUCTURE',
            'THERMAL_STORAGE': 'THERMAL STORAGE',
            'ELECTRIC_STORAGE': 'ELECTRIC STORAGE'
        }
        
        df_full['Category_Display'] = df_full['Category'].map(lambda x: category_display_map.get(x, x) if pd.notna(x) else x)
        
        # Count unmapped technologies
        unmapped = df_full[df_full['Category'].isna()]['Technologies'].unique()
        if len(unmapped) > 0:
            print(f"\nWarning: {len(unmapped)} technologies not mapped to categories:")
            for tech in unmapped[:10]:  # Show first 10
                print(f"  - {tech}")
            if len(unmapped) > 10:
                print(f"  ... and {len(unmapped) - 10} more")
        
        print(f"\nCategories found: {df_full['Category_Display'].nunique()}")
        
        # Create output directory
        output_dir = Path(output_dir)
        graphs_dir = output_dir / 'New_old_decom_graphs'
        graphs_dir.mkdir(parents=True, exist_ok=True)
        
        # Save full data CSV (like the example)
        full_csv = graphs_dir / 'New_old_decom_full_data.csv'
        df_full.to_csv(full_csv, index=False, sep=';')
        print(f"\n✓ Full data saved: {full_csv}")
        
        # Color mapping
        colors = px.colors.qualitative.Plotly + px.colors.qualitative.Set2 + px.colors.qualitative.Pastel
        unique_techs = sorted(df_full['Technologies'].unique())
        color_map = {tech: colors[i % len(colors)] for i, tech in enumerate(unique_techs)}
        
        # Generate plots by category
        print(f"\nGenerating plots by category...")
        categories_processed = 0
        
        for category in df_full['Category_Display'].dropna().unique():
            category_data = df_full[df_full['Category_Display'] == category].copy()
            
            if category_data.empty:
                continue
            
            category_name = str(category).replace('/', '_').replace(' ', '_').replace('&', 'and')
            category_dir = graphs_dir / category_name
            category_dir.mkdir(exist_ok=True)
            
            # Save category CSV
            category_csv = category_dir / f'{category_name}_data.csv'
            category_data.to_csv(category_csv, index=False, sep=';')
            
            # Phase order (historical first if present)
            phase_order = []
            if '1931_2015' in category_data['Phases'].unique():
                phase_order.append('1931_2015')
            phase_order.extend(sorted([p for p in category_data['Phases'].unique() 
                                       if p != '1931_2015']))
            
            try:
                if is_multiregional:
                    # Multi-regional plot
                    fig = px.bar(category_data,
                                x='Phases',
                                y='Tech_Cap',
                                color='Technologies',
                                facet_row='Type',
                                facet_col='Regions',
                                title=f'{category} - New/Old/Decom by Region',
                                color_discrete_map=color_map,
                                category_orders={'Phases': phase_order})
                    
                    fig.update_xaxes(tickangle=45)
                    fig.update_yaxes(matches=None, title_text="Capacity [GW]")
                    fig.update_layout(height=800, showlegend=True)
                    
                    html_path = category_dir / f'{category_name}_all_regions.html'
                    fig.write_html(str(html_path))
                    
                    try:
                        pdf_path = category_dir / f'{category_name}_all_regions.pdf'
                        fig.write_image(str(pdf_path), width=1600, height=800)
                    except:
                        pass
                    
                    # Individual region plots
                    for region in category_data['Regions'].unique():
                        region_data = category_data[category_data['Regions'] == region]
                        
                        fig_region = px.bar(region_data,
                                           x='Phases',
                                           y='Tech_Cap',
                                           color='Technologies',
                                           facet_row='Type',
                                           title=f'{category} - {region} - New/Old/Decom',
                                           color_discrete_map=color_map,
                                           category_orders={'Phases': phase_order})
                        
                        fig_region.update_xaxes(tickangle=45)
                        fig_region.update_yaxes(matches=None, title_text="Capacity [GW]")
                        fig_region.update_layout(height=600, showlegend=True)
                        
                        region_dir = category_dir / str(region)
                        region_dir.mkdir(exist_ok=True)
                        
                        html_path = region_dir / f'{category_name}_{region}.html'
                        fig_region.write_html(str(html_path))
                    
                else:
                    # Single-region plot
                    fig = px.bar(category_data,
                                x='Phases',
                                y='Tech_Cap',
                                color='Technologies',
                                facet_row='Type',
                                title=f'{category} - New/Old/Decom',
                                color_discrete_map=color_map,
                                category_orders={'Phases': phase_order})
                    
                    fig.update_xaxes(tickangle=45)
                    fig.update_yaxes(matches=None, title_text="Capacity [GW]")
                    fig.update_layout(height=600, showlegend=True)
                    
                    html_path = category_dir / f'{category_name}.html'
                    fig.write_html(str(html_path))
                    
                    try:
                        pdf_path = category_dir / f'{category_name}.pdf'
                        fig.write_image(str(pdf_path), width=1200, height=600)
                    except:
                        pass
                
                categories_processed += 1
                print(f"   {category} ({categories_processed}/{df_full['Category_Display'].nunique()})")
                
            except Exception as e:
                print(f"  ✗ Error plotting {category}: {e}")
        
        # Summary plot
        print("\nGenerating summary plot...")
        
        phase_order_all = []
        if '1931_2015' in df_full['Phases'].unique():
            phase_order_all.append('1931_2015')
        phase_order_all.extend(sorted([p for p in df_full['Phases'].unique() 
                                       if p != '1931_2015']))
        
        try:
            if is_multiregional:
                fig_summary = px.bar(df_full,
                                    x='Phases',
                                    y='Tech_Cap',
                                    color='Technologies',
                                    facet_row='Type',
                                    facet_col='Regions',
                                    title='All Categories - New/Old/Decom',
                                    color_discrete_map=color_map,
                                    category_orders={'Phases': phase_order_all})
                
                fig_summary.update_xaxes(tickangle=45)
                fig_summary.update_yaxes(matches=None, title_text="Capacity [GW]")
                fig_summary.update_layout(height=1000, showlegend=True)
                
                html_path = graphs_dir / 'summary_all_categories_regions.html'
                fig_summary.write_html(str(html_path))
            else:
                fig_summary = px.bar(df_full,
                                    x='Phases',
                                    y='Tech_Cap',
                                    color='Technologies',
                                    facet_row='Type',
                                    title='All Categories - New/Old/Decom',
                                    color_discrete_map=color_map,
                                    category_orders={'Phases': phase_order_all})
                
                fig_summary.update_xaxes(tickangle=45)
                fig_summary.update_yaxes(matches=None, title_text="Capacity [GW]")
                fig_summary.update_layout(height=800, showlegend=True)
                
                html_path = graphs_dir / 'summary_all_categories.html'
                fig_summary.write_html(str(html_path))
            
            print(f"✓ Summary saved")
        except Exception as e:
            print(f"Warning: Could not generate summary plot: {e}")
        
        print(f"\n{'='*70}")
        print(f"✓ Processing complete!")
        print(f"  Output directory: {graphs_dir}")
        print(f"  Categories processed: {categories_processed}")
        print(f"  Total data points: {len(df_full)}")
        print(f"{'='*70}\n")
        
        return df_full

    
        # ============================================================================
    # METHOD 1: Investment Cost Plotting
    # ============================================================================

    def graph_cost_inv_phase_tech(self, plot=True):
        """
        Plot cumulative investment costs by technology category over phases
        
        Aggregates C_inv_phase_tech over REGIONS and maps technologies to display categories.
        Creates cumulative area plots showing investment evolution across the pathway.
        
        Parameters
        ----------
        plot : bool, optional
            If True, generates and saves HTML/PDF plots (default: True)
            
        Returns
        -------
        pd.DataFrame
            Full plotting data with columns: Years, Category, C_inv_phase_tech, cumsum
        """
        logging.info("Generating investment cost plots by technology category")
        
        # Get C_inv_phase_tech from AMPL results
        # Dimensions: PHASE, REGIONS, TECHNOLOGIES
        try:
            results = self.esom.get_var('C_inv_phase_tech').copy()
        except Exception as e:
            logging.error(f"Could not retrieve C_inv_phase_tech: {e}")
            return None
        
        if results.empty:
            logging.warning("C_inv_phase_tech is empty, skipping plot")
            return None
        
        # Aggregate over REGIONS (sum across all regions)
        results_grouped = results.groupby(['Phase', 'Technologies']).sum()
        results_grouped.reset_index(inplace=True)
        
        # Filter out low values
        results_grouped = results_grouped[results_grouped['C_inv_phase_tech'] > 0.01]
        
        # Map technologies to display categories
        tech_category_map = self._get_tech_to_category_mapping()
        results_grouped['Category'] = results_grouped['Technologies'].map(tech_category_map)
        
        # Drop unmapped technologies
        unmapped = results_grouped[results_grouped['Category'].isna()]
        if len(unmapped) > 0:
            logging.warning(f"Dropping {len(unmapped)} unmapped technologies")
            logging.debug(f"Unmapped: {unmapped['Technologies'].unique().tolist()}")
        results_grouped = results_grouped.dropna(subset=['Category'])
        
        # Extract year from phase name (e.g., '2021_2025' -> '2025')
        results_grouped['Years'] = results_grouped['Phase'].map(lambda x: x.split('_')[-1])
        
        # Apply display name mapping
        category_display_map = {
            'ELECTRICITY': 'ELECTRICITY',
            'LIGHTING': 'LIGHTING',
            'HEAT_HIGH_T': 'HEAT HIGH TEMPERATURE',
            'HEAT_LOW_T': 'HEAT LOW TEMPERATURE',
            'COOLING': 'COOLING & REFRIGERATION',
            'FOOD_PRESERVATION': 'COOLING & REFRIGERATION',
            'COOKING': 'COOKING',
            'MOBILITY_PASSENGER': 'MOBILITY PASSENGER',
            'MOBILITY_FREIGHT': 'MOBILITY FREIGHT',
            'MECHANICAL_ENERGY': 'MECHANICAL ENERGY',
            'NON_ENERGY': 'NON-ENERGY',
            'BIOFUELS_EFUELS': 'BIOFUELS & EFUELS PRODUCTION',
            'ELEC_INFRASTRUCTURE': 'ELECTRICITY INFRASTRUCTURE',
            'HYDROCARBON_INFRASTRUCTURE': 'HYDROCARBON INFRASTRUCTURE',
            'H2_INFRASTRUCTURE': 'H2 INFRASTRUCTURE',
            'OTHER_INFRASTRUCTURE': 'OTHER INFRASTRUCTURE',
            'THERMAL_STORAGE': 'THERMAL STORAGE',
            'ELECTRIC_STORAGE': 'ELECTRIC STORAGE'
        }
        
        results_grouped['Category'] = results_grouped['Category'].map(lambda x: category_display_map.get(x, x))
        
        # Group by Years and Category
        df_to_plot = results_grouped.groupby(['Years', 'Category']).sum(numeric_only=True).reset_index()
        
        # Create full matrix (all year-category combinations with 0 for missing)
        years = sorted(df_to_plot['Years'].unique())
        categories = df_to_plot['Category'].unique()
        
        mi_temp = pd.MultiIndex.from_product([years, categories], names=['Years', 'Category'])
        temp = pd.DataFrame(0, index=mi_temp, columns=['C_inv_phase_tech'])
        temp.reset_index(inplace=True)
        
        df_to_plot = pd.concat([temp, df_to_plot])
        df_to_plot = df_to_plot.set_index(['Years', 'Category'])
        df_to_plot = df_to_plot.loc[~df_to_plot.index.duplicated(keep='last')]
        df_to_plot.sort_values(by=['Years'], inplace=True)
        
        # Calculate cumulative sum for each category
        for g_name, g_df in df_to_plot.groupby(['Category']):
            df_to_plot.loc[g_df.index, 'cumsum'] = df_to_plot.loc[g_df.index, 'C_inv_phase_tech'].cumsum()
        
        df_to_plot.reset_index(inplace=True)
        df_to_plot['cumsum'] = df_to_plot['cumsum'] / 1000  # Convert to billions
        df_to_plot['C_inv_phase_tech'] = df_to_plot['C_inv_phase_tech'] / 1000
        
        # Define category order (bottom to top for area chart)
        order_entry = [
            'OTHER INFRASTRUCTURE', 'H2 INFRASTRUCTURE', 'HYDROCARBON INFRASTRUCTURE',
            'ELECTRICITY INFRASTRUCTURE', 'BIOFUELS & EFUELS PRODUCTION',
            'THERMAL STORAGE', 'ELECTRIC STORAGE',
            'MOBILITY PASSENGER', 'MOBILITY FREIGHT',
            'ELECTRICITY', 'LIGHTING', 'HEAT HIGH TEMPERATURE', 'HEAT LOW TEMPERATURE',
            'COOLING & REFRIGERATION', 'NON-ENERGY', 'COOKING', 'MECHANICAL ENERGY'
        ]
        
        order_entry = [x for x in order_entry if x in list(df_to_plot['Category'])]
        
        df_to_plot['Category'] = pd.Categorical(df_to_plot['Category'], order_entry)
        df_to_plot.sort_values(by=['Category', 'Years'], inplace=True)
        
        # Define color dictionary
        cost_inv_color_dict = {
            "MOBILITY PASSENGER": "seagreen",
            "MOBILITY FREIGHT": "#8c564b",
            "ELECTRICITY": "#1f77b4",
            "HEAT HIGH TEMPERATURE": "#d62728",
            "COOKING": "#9467bd",
            "MECHANICAL ENERGY": "#7f7f7f",
            "NON-ENERGY": "darkgrey",
            "LIGHTING": "#bcbd22",
            "COOLING & REFRIGERATION": "#17becf",
            "HEAT LOW TEMPERATURE": "fuchsia",
            "BIOFUELS & EFUELS PRODUCTION": "orange",
            "ELECTRICITY INFRASTRUCTURE": "lightblue",
            "HYDROCARBON INFRASTRUCTURE": "brown",
            "H2 INFRASTRUCTURE": "cyan",
            "OTHER INFRASTRUCTURE": "grey",
            "THERMAL STORAGE": "pink",
            "ELECTRIC STORAGE": "chartreuse"
        }
        
        # Save CSV file
        csv_file_path = Path(self.cs_dir) / "outputs" / "Cost_Inv_Phase" / "cost_inv_phase_tech.csv"
        csv_file_path.parent.mkdir(parents=True, exist_ok=True)
        df_to_plot.to_csv(csv_file_path, index=False)
        logging.info(f"Saved CSV: {csv_file_path}")
        
        if plot:
            try:
                import plotly.express as px
            except ImportError:
                logging.warning("Plotly not installed. Skipping plot generation, but CSV is saved")
                return df_to_plot
            
            # Create area plot
            fig = px.area(
                df_to_plot, 
                x='Years', 
                y='cumsum', 
                color='Category',
                title=self.case_study + ' - Investment Costs',
                color_discrete_map=cost_inv_color_dict,
                category_orders={"Category": order_entry}
            )
            
            fig.for_each_trace(lambda trace: trace.update(fillcolor=trace.line.color))
            fig.update_traces(mode='none')
            fig.update_xaxes(categoryorder='array', categoryarray=sorted(df_to_plot['Years'].unique()))
            
            # Layout adjustments
            title = "<b>" + self.case_study + "</b><br>[Investments over transition - bâ‚¬<sub>2021</sub>]"
            fig.update_layout(
                title=dict(text=title, x=0.5),
                yaxis_title='Cumulative Investment (bâ‚¬<sub>2021</sub>)',
                legend=dict(
                    orientation="v",
                    yanchor="top",
                    y=1,
                    xanchor="left",
                    x=1.04,
                    traceorder="reversed"
                )
            )
            
            # Create output directory
            outdir_path = Path(self.cs_dir) / "outputs" / "Cost_Inv_Phase"
            outdir_path.mkdir(parents=True, exist_ok=True)
            
            # Save as HTML
            html_file_path = str(outdir_path / "cost_inv_phase_tech.html")
            fig.write_html(html_file_path)
            
            # Save as PDF
            try:
                pdf_file_path = str(outdir_path / "cost_inv_phase_tech.pdf")
                fig.write_image(pdf_file_path, width=1200, height=550)
            except Exception as e:
                logging.warning(f"Could not save PDF: {e}")
            
            # Also save to main output directory
            html_file_path_main = str(self.cs_dir / "outputs" / "final_cost_inv_phase_tech.html")
            fig.write_html(html_file_path_main)
            try:
                pdf_file_path_main = str(self.cs_dir / "outputs" / "final_cost_inv_phase_tech.pdf")
                fig.write_image(pdf_file_path_main, width=1200, height=550)
            except:
                pass
            
            logging.info(f"Saved plots: {outdir_path}")
        
        return df_to_plot


    # ============================================================================
    # METHOD 2: Operation Cost Plotting
    # ============================================================================

    def graph_cost_op_phase(self, plot=True):
        """
        Plot cumulative operation costs by category over phases
        
        Combines C_op_phase_tech (technologies) and C_op_phase_res (resources),
        aggregates over REGIONS, separates RE vs NRE fuels, and creates cumulative plots.
        
        Parameters
        ----------
        plot : bool, optional
            If True, generates and saves HTML/PDF plots (default: True)
            
        Returns
        -------
        pd.DataFrame
            Full plotting data with columns: Years, Category, C_op, cumsum
        """
        logging.info("Generating operation cost plots by category")
        
        # Get C_op_phase_tech (technologies)
        try:
            results_tech = self.esom.get_var('C_op_phase_tech').copy()
        except Exception as e:
            logging.error(f"Could not retrieve C_op_phase_tech: {e}")
            return None
        
        # Get C_op_phase_res (resources)
        try:
            results_res = self.esom.get_var('C_op_phase_res').copy()
        except Exception as e:
            logging.error(f"Could not retrieve C_op_phase_res: {e}")
            results_res = pd.DataFrame()
        
        if results_tech.empty:
            logging.warning("C_op_phase_tech is empty, skipping plot")
            return None
        
        # Aggregate over REGIONS
        results_tech_grouped = results_tech.groupby(['Phase', 'Technologies']).sum()
        results_tech_grouped.reset_index(inplace=True)
        results_tech_grouped.rename(columns={'C_op_phase_tech': 'C_op', 'Technologies': 'Elements'}, inplace=True)
        results_tech_grouped['Years'] = results_tech_grouped['Phase'].map(lambda x: x.split('_')[-1])
        
        if not results_res.empty:
            results_res_grouped = results_res.groupby(['Phase', 'Resources']).sum()
            results_res_grouped.reset_index(inplace=True)
            results_res_grouped.rename(columns={'C_op_phase_res': 'C_op', 'Resources': 'Elements'}, inplace=True)
            results_res_grouped['Years'] = results_res_grouped['Phase'].map(lambda x: x.split('_')[-1])
        else:
            results_res_grouped = pd.DataFrame(columns=['Phase', 'Elements', 'C_op', 'Years'])
        
        # Determine first operational year to avoid duplication with Cost_breakdown
        if len(self.years) > 1:
            first_op_year = str(self.years[1])  # e.g., '2021' in a 2015-2050 pathway
        else:
            first_op_year = None
        
        # Filter out first operational year from phase data (will be added from Cost_breakdown)
        if first_op_year:
            results_tech_grouped = results_tech_grouped[results_tech_grouped['Years'] != first_op_year]
            results_res_grouped = results_res_grouped[results_res_grouped['Years'] != first_op_year]
            logging.info(f"Filtered out {first_op_year} from phase data to avoid duplication with Cost_breakdown")
        
        # Get Cost_breakdown for first operational year (e.g., YEAR_2021)
        # Note: In pathway models, years[0] is base year (2015), years[1] is first operational year (2021)
        if first_op_year:
            first_op_year_key = f'YEAR_{first_op_year}'
            
            try:
                cost_breakdown = self.results['Cost_breakdown'].copy()
                # Filter for first operational year and aggregate over regions
                results_first_year = cost_breakdown.loc[
                    cost_breakdown.index.get_level_values('Years') == first_op_year_key, :
                ]
                
                if not results_first_year.empty:
                    results_first_year = results_first_year.groupby('Elements').sum()
                    results_first_year.reset_index(inplace=True)
                    
                    # Combine C_op and C_maint for first year
                    if 'C_op' in results_first_year.columns and 'C_maint' in results_first_year.columns:
                        results_first_year['C_op'] = results_first_year['C_op'].fillna(0) + results_first_year['C_maint'].fillna(0)
                    elif 'C_maint' in results_first_year.columns:
                        results_first_year['C_op'] = results_first_year['C_maint'].fillna(0)
                    
                    # Keep only Elements and C_op columns
                    results_first_year = results_first_year[['Elements', 'C_op']]
                    results_first_year['Years'] = first_op_year
                    
                    logging.info(f"Added Cost_breakdown data for first operational year {first_op_year}")
                else:
                    results_first_year = pd.DataFrame(columns=['Elements', 'C_op', 'Years'])
                    logging.warning(f"No Cost_breakdown data found for {first_op_year_key}")
            except Exception as e:
                logging.warning(f"Could not retrieve Cost_breakdown for first operational year: {e}")
                results_first_year = pd.DataFrame(columns=['Elements', 'C_op', 'Years'])
        else:
            results_first_year = pd.DataFrame(columns=['Elements', 'C_op', 'Years'])
        
        # Combine first year, tech and resource costs
        df_to_plot = pd.concat([results_first_year, results_tech_grouped, results_res_grouped], 
                               axis=0, ignore_index=True)
        
        # Map elements to categories
        tech_category_map = self._get_tech_to_category_mapping()
        df_to_plot['Category'] = df_to_plot['Elements'].map(tech_category_map)
        
        # For unmapped items, try resource mapping
        unmapped_mask = df_to_plot['Category'].isna()
        if unmapped_mask.any():
            resource_category_map = self._get_resource_to_category_mapping()
            df_to_plot.loc[unmapped_mask, 'Category'] = df_to_plot.loc[unmapped_mask, 'Elements'].map(resource_category_map)
        
        # Drop remaining unmapped
        unmapped = df_to_plot[df_to_plot['Category'].isna()]
        if len(unmapped) > 0:
            logging.warning(f"Dropping {len(unmapped)} unmapped elements")
            logging.debug(f"Unmapped: {unmapped['Elements'].unique().tolist()}")
            
            # Save unmapped elements to CSV
            unmapped_csv_path = Path(self.cs_dir) / "outputs" / "Cost_Op_Phase" / "unmapped_elements.csv"
            unmapped_csv_path.parent.mkdir(parents=True, exist_ok=True)
            
            # Create DataFrame with unmapped element details
            unmapped_summary = unmapped[['Elements', 'Years', 'C_op']].copy()
            unmapped_summary = unmapped_summary.groupby('Elements').agg({
                'Years': lambda x: ', '.join(sorted(set(x))),
                'C_op': 'sum'
            }).reset_index()
            unmapped_summary.columns = ['Element', 'Years_Present', 'Total_C_op']
            unmapped_summary = unmapped_summary.sort_values('Total_C_op', ascending=False)
            unmapped_summary.to_csv(unmapped_csv_path, index=False)
            
            logging.info(f"Saved unmapped elements list to: {unmapped_csv_path}")
        
        df_to_plot = df_to_plot.dropna(subset=['Category'])
        
        # Apply display name mapping
        category_display_map = {
            'ELECTRICITY': 'ELECTRICITY',
            'LIGHTING': 'LIGHTING',
            'HEAT_HIGH_T': 'HEAT HIGH TEMPERATURE',
            'HEAT_LOW_T': 'HEAT LOW TEMPERATURE',
            'COOLING': 'COOLING & REFRIGERATION',
            'FOOD_PRESERVATION': 'COOLING & REFRIGERATION',
            'COOKING': 'COOKING',
            'MOBILITY_PASSENGER': 'MOBILITY PASSENGER',
            'MOBILITY_FREIGHT': 'MOBILITY FREIGHT',
            'MECHANICAL_ENERGY': 'MECHANICAL ENERGY',
            'NON_ENERGY': 'NON-ENERGY',
            'BIOFUELS_EFUELS': 'BIOFUELS & EFUELS PRODUCTION',
            'ELEC_INFRASTRUCTURE': 'ELECTRICITY INFRASTRUCTURE',
            'HYDROCARBON_INFRASTRUCTURE': 'HYDROCARBON INFRASTRUCTURE',
            'H2_INFRASTRUCTURE': 'H2 INFRASTRUCTURE',
            'OTHER_INFRASTRUCTURE': 'OTHER INFRASTRUCTURE',
            'THERMAL_STORAGE': 'THERMAL STORAGE',
            'ELECTRIC_STORAGE': 'ELECTRIC STORAGE',
            'RE_FUELS': 'RE_FUELS',
            'NRE_FUELS': 'NRE_FUELS'
        }
        
        df_to_plot['Category'] = df_to_plot['Category'].map(lambda x: category_display_map.get(x, x))
        
        # Group by Years and Category
        df_to_plot = df_to_plot.groupby(['Years', 'Category']).sum(numeric_only=True).reset_index()
        
        # Create full matrix
        years = sorted(df_to_plot['Years'].unique())
        categories = df_to_plot['Category'].unique()
        
        mi_temp = pd.MultiIndex.from_product([years, categories], names=['Years', 'Category'])
        temp = pd.DataFrame(0, index=mi_temp, columns=['C_op'])
        temp.reset_index(inplace=True)
        
        df_to_plot = pd.concat([temp, df_to_plot])
        df_to_plot = df_to_plot.set_index(['Years', 'Category'])
        df_to_plot = df_to_plot.loc[~df_to_plot.index.duplicated(keep='last')]
        df_to_plot.sort_values(by=['Years'], inplace=True)
        
        # Calculate cumulative sum
        for g_name, g_df in df_to_plot.groupby(['Category']):
            df_to_plot.loc[g_df.index, 'cumsum'] = df_to_plot.loc[g_df.index, 'C_op'].cumsum()
        
        df_to_plot.reset_index(inplace=True)
        df_to_plot['cumsum'] = df_to_plot['cumsum'] / 1000  # Convert to billions
        df_to_plot['C_op'] = df_to_plot['C_op'] / 1000
        
        # Define category order (with fuel types at bottom)
        order_entry = [
            'NRE_FUELS', 'RE_FUELS',
            'OTHER INFRASTRUCTURE', 'H2 INFRASTRUCTURE', 'HYDROCARBON INFRASTRUCTURE',
            'ELECTRICITY INFRASTRUCTURE', 'BIOFUELS & EFUELS PRODUCTION',
            'THERMAL STORAGE', 'ELECTRIC STORAGE',
            'MOBILITY PASSENGER', 'MOBILITY FREIGHT',
            'ELECTRICITY', 'LIGHTING', 'HEAT HIGH TEMPERATURE', 'HEAT LOW TEMPERATURE',
            'COOLING & REFRIGERATION', 'NON-ENERGY', 'COOKING', 'MECHANICAL ENERGY'
        ]
        
        order_entry = [x for x in order_entry if x in list(df_to_plot['Category'])]
        
        df_to_plot['Category'] = pd.Categorical(df_to_plot['Category'], order_entry)
        df_to_plot.sort_values(by=['Category', 'Years'], inplace=True)
        
        # Define color dictionary
        cost_op_color_dict = {
            "MOBILITY PASSENGER": "seagreen",
            "MOBILITY FREIGHT": "#8c564b",
            "ELECTRICITY": "#1f77b4",
            "HEAT HIGH TEMPERATURE": "#d62728",
            "COOKING": "#9467bd",
            "MECHANICAL ENERGY": "#7f7f7f",
            "NON-ENERGY": "darkgrey",
            "LIGHTING": "#bcbd22",
            "COOLING & REFRIGERATION": "#17becf",
            "HEAT LOW TEMPERATURE": "fuchsia",
            "BIOFUELS & EFUELS PRODUCTION": "orange",
            "ELECTRICITY INFRASTRUCTURE": "lightblue",
            "HYDROCARBON INFRASTRUCTURE": "brown",
            "H2 INFRASTRUCTURE": "cyan",
            "OTHER INFRASTRUCTURE": "grey",
            "THERMAL STORAGE": "pink",
            "ELECTRIC STORAGE": "chartreuse",
            "RE_FUELS": "green",
            "NRE_FUELS": "black"
        }
        
        # Save CSV file
        csv_file_path = Path(self.cs_dir) / "outputs" / "Cost_Op_Phase" / "cost_op_phase.csv"
        csv_file_path.parent.mkdir(parents=True, exist_ok=True)
        df_to_plot.to_csv(csv_file_path, index=False)
        logging.info(f"Saved CSV: {csv_file_path}")
        
        if plot:
            try:
                import plotly.express as px
            except ImportError:
                logging.warning("Plotly not installed. Skipping plot generation, but CSV is saved")
                return df_to_plot
            
            # Create area plot
            fig = px.area(
                df_to_plot,
                x='Years',
                y='cumsum',
                color='Category',
                title=self.case_study + ' - Operation Costs',
                color_discrete_map=cost_op_color_dict,
                category_orders={"Category": order_entry}
            )
            
            fig.for_each_trace(lambda trace: trace.update(fillcolor=trace.line.color))
            fig.update_traces(mode='none')
            fig.update_xaxes(categoryorder='array', categoryarray=sorted(df_to_plot['Years'].unique()))
            
            # Layout adjustments
            title = "<b>" + self.case_study + "</b><br>[Operation over transition - bâ‚¬<sub>2021</sub>]"
            fig.update_layout(
                title=dict(text=title, x=0.5),
                yaxis_title='Cumulative Operation Cost (bâ‚¬<sub>2021</sub>)',
                legend=dict(
                    orientation="v",
                    yanchor="top",
                    y=1,
                    xanchor="left",
                    x=1.04,
                    traceorder="reversed"
                )
            )
            
            # Create output directory
            outdir_path = Path(self.cs_dir) / "outputs" / "Cost_Op_Phase"
            outdir_path.mkdir(parents=True, exist_ok=True)
            
            # Save as HTML
            html_file_path = str(outdir_path / "cost_op_phase.html")
            fig.write_html(html_file_path)
            
            # Save as PDF
            try:
                pdf_file_path = str(outdir_path / "cost_op_phase.pdf")
                fig.write_image(pdf_file_path, width=1200, height=550)
            except Exception as e:
                logging.warning(f"Could not save PDF: {e}")
            
            # Also save to main output directory
            html_file_path_main = str(self.cs_dir / "outputs" / "final_cost_op_phase.html")
            fig.write_html(html_file_path_main)
            try:
                pdf_file_path_main = str(self.cs_dir / "outputs" / "final_cost_op_phase.pdf")
                fig.write_image(pdf_file_path_main, width=1200, height=550)
            except:
                pass
            
            logging.info(f"Saved plots: {outdir_path}")
        
        return df_to_plot


    # ============================================================================
    # METHOD 3: Technology to Category Mapping
    # ============================================================================

    def _get_tech_to_category_mapping(self):
        """
        Create mapping from technology names to display categories
        
        Returns
        -------
        dict
            Mapping of technology name to category
        """
        # Based on END_USES_TYPES_OF_CATEGORY from the model
        tech_by_category = {
            "ELECTRICITY": [
                "NUCLEAR", "NUCLEAR_SMR", "CCGT_AMMONIA", "COAL_US", "COAL_IGCC", 
                "BIOMASS_TO_POWER", "PV_ROOFTOP", "PV_UTILITY", "PT_POWER_BLOCK", 
                "ST_POWER_BLOCK", "WIND_ONSHORE", "WIND_OFFSHORE", "HYDRO_DAM", 
                "HYDRO_RIVER", "TIDAL_STREAM", "TIDAL_RANGE", "WAVE", "GEOTHERMAL", 
                "GENSET_DIESEL", "OCGT", "CCGT", "ST_BIOMASS", "ST_SNG", "CCGT_SUR", 
                "FUEL_CELL", "FB_ST_BIOMASS", "CFB_ST_BIOMASS", "BFB_ST_BIOMASS"
            ],
            "LIGHTING": [
                "CONVENTIONAL_BULB", "LED_BULB", "CONVENTIONAL_LIGHT", "LED_LIGHT"
            ],
            "HEAT_HIGH_T": [
                "IND_COGEN_GAS", "IND_COGEN_WOOD", "IND_COGEN_WASTE", "IND_BOILER_GAS",
                "IND_BOILER_WOOD", "IND_BOILER_BIOWASTE", "IND_BOILER_OIL", "IND_BOILER_COAL",
                "IND_BOILER_WASTE", "IND_DIRECT_ELEC", "IND_BOILER_DIESEL", "IND_BOILER_LPG"
            ],
            "HEAT_LOW_T": [
                "DHN_HP_ELEC", "DHN_COGEN_GAS", "DHN_COGEN_WOOD", "DHN_COGEN_WASTE",
                "DHN_BOILER_GAS", "DHN_BOILER_WOOD", "DHN_BOILER_OIL", "DHN_DEEP_GEO",
                "DHN_SOLAR", "DEC_HP_ELEC", "DEC_THHP_GAS", "DEC_COGEN_GAS", "DEC_COGEN_OIL",
                "DEC_ADVCOGEN_GAS", "DEC_ADVCOGEN_H2", "DEC_BOILER_GAS", "DEC_BOILER_WOOD",
                "DEC_BOILER_OIL", "DEC_SOLAR", "DEC_DIRECT_ELEC"
            ],
            "COOLING": [
                "DEC_THHP_GAS_COLD", "DEC_ELEC_COLD", "IND_ELEC_COLD"
            ],
            "FOOD_PRESERVATION": [
                "REFRIGERATOR_EL"
            ],
            "COOKING": [
                "STOVE_WOOD", "STOVE_LPG", "STOVE_NG", "STOVE_OIL", "STOVE_ELEC"
            ],
            "MOBILITY_PASSENGER": [
                # MOB_PUBLIC
                "TRAMWAY_TROLLEY", "BUS_COACH_DIESEL", "BUS_COACH_HYDIESEL",
                "BUS_COACH_CNG_STOICH", "BUS_COACH_FC_HYBRIDH2", "TRAIN_PUB",
                "CAR_FG_PUBLIC", "BUS_FG_PUBLIC", "PICKUP_TRUCK_FG_PUBLIC",
                "SUV_FG_PUBLIC", "MOTORCYCLE_FG_PUBLIC", "CAR_GASOLINE_PUBLIC",
                "BUS_GASOLINE_PUBLIC", "PICKUP_TRUCK_GASOLINE_PUBLIC", "SUV_GASOLINE_PUBLIC",
                "MOTORCYCLE_GASOLINE_PUBLIC", "CAR_DIESEL_PUBLIC", "BUS_DIESEL_PUBLIC",
                "PICKUP_TRUCK_DIESEL_PUBLIC", "SUV_DIESEL_PUBLIC", "MOTORCYCLE_DIESEL_PUBLIC",
                "BUS_ELEC_PUBLIC", "PICKUP_TRUCK_ELEC_PUBLIC", "SUV_ELEC_PUBLIC",
                "MOTORCYCLE_ELEC_PUBLIC",
                # MOB_PRIVATE
                "CAR_GASOLINE", "CAR_DIESEL", "CAR_NG", "CAR_METHANOL", "CAR_HEV",
                "CAR_PHEV", "CAR_BEV", "CAR_FUEL_CELL", "CAR_FG_PRIVATE",
                "BUS_FG_PRIVATE", "PICKUP_TRUCK_FG_PRIVATE", "SUV_FG_PRIVATE",
                "MOTORCYCLE_FG_PRIVATE", "CAR_GASOLINE_PRIVATE", "BUS_GASOLINE_PRIVATE",
                "PICKUP_TRUCK_GASOLINE_PRIVATE", "SUV_GASOLINE_PRIVATE", "MOTORCYCLE_GASOLINE_PRIVATE",
                "CAR_DIESEL_PRIVATE", "BUS_DIESEL_PRIVATE", "PICKUP_TRUCK_DIESEL_PRIVATE",
                "SUV_DIESEL_PRIVATE", "MOTORCYCLE_DIESEL_PRIVATE", "BUS_ELEC_PRIVATE",
                "PICKUP_TRUCK_ELEC_PRIVATE", "SUV_ELEC_PRIVATE", "MOTORCYCLE_ELEC_PRIVATE"
            ],
            "MOBILITY_FREIGHT": [
                # MOB_FREIGHT_RAIL
                "TRAIN_FREIGHT", "TRAIN_FREIGHT_ELEC",
                # MOB_FREIGHT_BOAT
                "BOAT_FREIGHT_DIESEL", "BOAT_FREIGHT_NG", "BOAT_FREIGHT_METHANOL", "BOAT_FREIGHT_ELEC",
                # MOB_FREIGHT_ROAD
                "TRUCK_DIESEL", "TRUCK_METHANOL", "TRUCK_FUEL_CELL", "TRUCK_ELEC",
                "TRUCK_NG", "TRUCK_FG", "TRUCK_GASOLINE", "TRUCK_DIESEL_P",
                # AVIATION
                "PLANE_SHORT_HAUL", "PLANE_H2_SHORT_HAUL", "PLANE", "PLANE_LONG_HAUL",
                # SHIPPING
                "CARGO_LFO", "CARGO_LNG", "CARGO_METHANOL", "CARGO_AMMONIA",
                "CARGO_FUELCELL_LH2", "CARGO_FUELCELL_AMMONIA", "CARGO_RETRO_METHANOL",
                "CARGO_RETRO_AMMONIA", "VAN_FG", "CARGO_MOTORCYCLE_FG", "VAN_GASOLINE",
                "CARGO_MOTORCYCLE_GASOLINE", "VAN_DIESEL", "CARGO_MOTORCYCLE_DIESEL",
                "VAN_ELEC", "CARGO_MOTORCYCLE_ELEC"
            ],
            "MECHANICAL_ENERGY": [
                "COMM_MACHINERY_DIESEL", "COMM_MACHINERY_EL", "IND_MACHINERY_EL",
                "TRACTOR_DIESEL", "TRACTOR_EL", "AGR_MACHINERY_DIESEL", "AGR_MACHINERY_EL",
                "MIN_MACHINERY_DIESEL", "MIN_MACHINERY_EL", "FISH_MACHINERY_DIESEL",
                "FISH_MACHINERY_EL"
            ],
            "NON_ENERGY": [
                "HABER_BOSCH", "SYN_METHANOLATION", "METHANE_TO_METHANOL",
                "BIOMASS_TO_METHANOL", "BIOWASTE_TO_METHANOL", "OIL_TO_HVC",
                "GAS_TO_HVC", "BIOMASS_TO_HVC", "METHANOL_TO_HVC"
            ],
            "BIOFUELS_EFUELS": [
                "H2_ELECTROLYSIS", "H2_NG", "H2_BIOMASS", "BIOMASS_TO_METHANE",
                "BIOWASTE_TO_METHANE", "SYN_METHANATION", "BIOMETHANATION_WET_BIOMASS",
                "BIOMETHANATION_BIOWASTE", "BIOMASS_TO_GASOLINE", "BIOMASS_TO_DIESEL",
                "BIOMASS_TO_JET_FUEL", "BIOMASS_TO_LFO", "BIOWASTE_TO_GASOLINE",
                "BIOWASTE_TO_DIESEL", "BIOWASTE_TO_JET_FUEL", "BIOWASTE_TO_LFO",
                "DIESEL_TO_JET_FUEL", "ATM_CCS", "INDUSTRY_CCS", "POWER_TO_GASOLINE",
                "POWER_TO_DIESEL", "POWER_TO_JET_FUEL", "POWER_TO_LFO", "H2_TO_GASOLINE",
                "H2_TO_DIESEL", "H2_TO_JET_FUEL", "H2_TO_LFO", "AMMONIA_TO_H2",
                "FERMENTATION_TO_BIOETHANOL", "ESTERIFICATION_TO_BIODIESEL", "ETHANOL_TO_FUELS",
                "BIO_HYDROLYSIS", "PYROLYSIS_TO_FUELS", "PYROLYSIS_TO_LFO"
            ],
            "ELEC_INFRASTRUCTURE": [
                "HVAC_LINE", "HVDC_SUBSEA"
            ],
            "HYDROCARBON_INFRASTRUCTURE": [
                "GAS_PIPELINE", "GAS_SUBSEA", "DIESEL_PIPELINE", "GASOLINE_PIPELINE",
                "LPG_PIPELINE", "LFO_PIPELINE", "JET_FUEL_PIPELINE", "REGASIFICATION",
                "LIQUEFACTION"
            ],
            "H2_INFRASTRUCTURE": [
                "H2_RETROFITTED", "H2_NEW", "H2_SUBSEA_RETRO", "H2_SUBSEA_NEW"
            ],
            "OTHER_INFRASTRUCTURE": [
                "PT_COLLECTOR", "ST_COLLECTOR", "EFFICIENCY", "DHN", "GRID"
            ],
            "THERMAL_STORAGE": [
                "TS_DEC_HP_ELEC", "TS_DEC_THHP_GAS", "TS_DEC_COGEN_GAS", "TS_DEC_COGEN_OIL",
                "TS_DEC_ADVCOGEN_GAS", "TS_DEC_ADVCOGEN_H2", "TS_DEC_BOILER_GAS",
                "TS_DEC_BOILER_WOOD", "TS_DEC_BOILER_OIL", "TS_DEC_DIRECT_ELEC",
                "TS_DHN_DAILY", "TS_DHN_SEASONAL", "TS_HIGH_TEMP", "TS_COLD"
            ],
            "ELECTRIC_STORAGE": [
                "DAM_STORAGE", "PHS", "BATT_LI", "BEV_BATT", "PHEV_BATT", "CAES",
                "SUV_ELEC_PRIVATE_BATT", "MOTORCYCLE_ELEC_PRIVATE_BATT", "BUS_ELEC_PUBLIC_BATT"
            ]
        }
        
        # Create reverse mapping (tech -> category)
        tech_to_cat = {}
        for cat, techs in tech_by_category.items():
            for tech in techs:
                tech_to_cat[tech] = cat
        
        # Add storage technologies that might be mapped differently
        for tech in ["GAS_STORAGE", "H2_STORAGE", "DIESEL_STORAGE", "JET_FUEL_STORAGE",
                    "GASOLINE_STORAGE", "LFO_STORAGE", "AMMONIA_STORAGE", "METHANOL_STORAGE",
                    "CO2_STORAGE", "PT_STORAGE", "ST_STORAGE", "LPG_STORAGE", "LNG_STORAGE"]:
            if tech not in tech_to_cat:
                if tech in ["PT_STORAGE", "ST_STORAGE", "TS_COLD"]:
                    tech_to_cat[tech] = "THERMAL_STORAGE"
                else:
                    tech_to_cat[tech] = "OTHER_INFRASTRUCTURE"
        
        return tech_to_cat


    # ============================================================================
    # METHOD 4: Resource to Category Mapping
    # ============================================================================

    def _get_resource_to_category_mapping(self):
        """
        Create mapping from resource names to display categories
        
        For operation costs, resources need to be categorized as RE_FUELS or NRE_FUELS,
        or mapped to other appropriate categories.
        
        Returns
        -------
        dict
            Mapping of resource name to category
        """
        # RE fuels
        re_fuels = [
            "GASOLINE_RE", "DIESEL_RE", "LFO_RE", "JET_FUEL_RE", "GAS_RE",
            "H2_RE", "AMMONIA_RE", "METHANOL_RE", "WOOD", "WET_BIOMASS",
            "ENERGY_CROPS_2", "BIOMASS_RESIDUES", "BIOWASTE"
        ]
        
        # NRE fuels
        nre_fuels = [
            "GASOLINE", "DIESEL", "LFO", "JET_FUEL", "GAS", "H2", "AMMONIA",
            "METHANOL", "COAL", "URANIUM", "WASTE", "LPG", "LNG"
        ]
        
        # Electricity is its own category
        electricity = ["ELECTRICITY"]
        
        resource_to_cat = {}
        
        for res in re_fuels:
            resource_to_cat[res] = "RE_FUELS"
        
        for res in nre_fuels:
            resource_to_cat[res] = "NRE_FUELS"
        
        for res in electricity:
            resource_to_cat[res] = "ELECTRICITY"
        
        return resource_to_cat