#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate Individual Year .dat Files
====================================

This script generates the seven .dat files for each year from the existing
data structure, using the same methodology as the one-year model.

Files generated for each year:
  1. indep.dat
  2. reg_demands.dat
  3. reg_resources.dat
  4. reg_technologies.dat
  5. reg_exch.dat
  6. reg_misc.dat
  7. reg_storage_power_to_energy.dat

Usage:
    python generate_year_dat_files.py

The script will process all year directories found in the Data directory
(2015, 2021, 2025, 2030, 2035, 2040, 2045, 2050) and generate the .dat
files in each respective year directory.
"""

import logging
import copy
import numpy as np
import pandas as pd
import sys
from pathlib import Path

# Add project path
sys.path.append('/home/pjimenez/ESMC_PTH/EnergyScope_Unif/')

from esmc.utils.region import Region
from esmc.utils.df_utils import clean_indices
import esmc.preprocessing.dat_print as dp
from esmc.common import CSV_SEPARATOR, bo_country_code, named_space_id
import esmc.postprocessing.amplpy2pd as a2p

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%H:%M:%S'
)


class YearDataGenerator:
    """
    Generator for individual year .dat files
    Uses the same logic as the one-year Esmc model
    """
    
    def __init__(self, year, regions_names, data_dir, output_dir=None):
        """
        Initialize the year data generator
        
        Parameters
        ----------
        year : int or str
            Year to process
        regions_names : list
            List of region names
        data_dir : Path
            Base data directory
        output_dir : Path, optional
            Output directory for .dat files. If None, saves in year directory.
        """
        self.year = str(year)
        self.regions_names = regions_names
        self.regions_names.sort()
        
        # Paths
        self.project_dir = Path(__file__).parents[0]
        self.data_dir = data_dir
        self.year_dir = data_dir / self.year
        
        if output_dir is None:
            self.output_dir = self.year_dir
        else:
            self.output_dir = output_dir
            
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Data containers
        self.ref_region_name = '02_REF_REGION'
        self.ref_region = None
        self.regions = dict.fromkeys(self.regions_names, None)
        
        self.sets = dict()
        self.data_indep = dict.fromkeys([
            'Misc_indep', 'Layers_in_out', 'Resources_indep',
            'Storage_characteristics', 'Storage_eff_in', 'Storage_eff_out'
        ])
        self.data_reg = dict.fromkeys([
            'Exch', 'Demands', 'Technologies',
            'Storage_power_to_energy', 'Misc'
        ])
        
        logging.info(f'Initialized YearDataGenerator for year {self.year}')
    
    def read_data_indep(self):
        """Read data independent of the region dimension"""
        data_path = self.year_dir / '00_INDEP'
        
        if not data_path.exists():
            raise FileNotFoundError(f'INDEP data directory not found: {data_path}')
        
        logging.info(f'Reading indep data from {data_path}')
        
        # Reading misc_indep
        r_path = data_path / 'Misc_indep.json'
        if r_path.is_file():
            self.data_indep['Misc_indep'] = a2p.read_json(r_path)
        else:
            logging.warning(f'Misc_indep.json not found at {r_path}')
            self.data_indep['Misc_indep'] = {}
        
        # Reading layers_in_out
        self.data_indep['Layers_in_out'] = pd.read_csv(
            data_path / 'Layers_in_out.csv',
            sep=CSV_SEPARATOR,
            header=[0],
            index_col=[0]
        )
        self.data_indep['Layers_in_out'] = clean_indices(self.data_indep['Layers_in_out'])
        
        # Reading resources_indep
        self.data_indep['Resources_indep'] = clean_indices(
            pd.read_csv(
                data_path / 'Resources_indep.csv',
                sep=CSV_SEPARATOR,
                header=[2],
                index_col=[2]
            )
        ).drop(columns=['Comment'], errors='ignore')
        
        # Reading storage_characteristics
        self.data_indep['Storage_characteristics'] = clean_indices(
            pd.read_csv(
                data_path / 'Storage_characteristics.csv',
                sep=CSV_SEPARATOR,
                header=[0],
                index_col=[0]
            )
        )
        
        # Reading storage_eff_in
        self.data_indep['Storage_eff_in'] = clean_indices(
            pd.read_csv(
                data_path / 'Storage_eff_in.csv',
                sep=CSV_SEPARATOR,
                header=[0],
                index_col=[0]
            )
        ).dropna(axis=0, how='all')
        
        # Reading storage_eff_out
        self.data_indep['Storage_eff_out'] = clean_indices(
            pd.read_csv(
                data_path / 'Storage_eff_out.csv',
                sep=CSV_SEPARATOR,
                header=[0],
                index_col=[0]
            )
        ).dropna(axis=0, how='all')
        
        return
    
    def init_regions(self):
        """Initialize regions and read their data"""
        logging.info(f'Initializing regions: {", ".join(self.regions_names)}')
        
        # Initialize reference region
        self.ref_region = Region(
            nuts=self.ref_region_name,
            data_dir=self.year_dir,
            ref_region=True
        )
        
        # Initialize other regions
        for r in self.regions_names:
            if r != self.ref_region_name:
                self.regions[r] = copy.deepcopy(self.ref_region)
                self.regions[r].__init__(nuts=r, data_dir=self.year_dir, ref_region=False)
            else:
                self.regions[r] = self.ref_region
        
        # Read exchange data
        self.read_data_exch()
        
        return
    
    def read_data_exch(self):
        """Read exchange data between regions"""
        data_path = self.year_dir
        
        # Initialize Exch data structure
        self.data_reg['Exch'] = {}
        
        # Read Misc_exch
        r_path = data_path / '01_EXCH' / 'Misc_exch.json'
        if r_path.is_file():
            self.data_reg['Exch']['Misc_exch'] = a2p.read_json(r_path)
        else:
            logging.warning(f'Misc_exch.json not found, using empty dict')
            self.data_reg['Exch']['Misc_exch'] = {'add_sets': {}}
        
        # Try to read exchange network data if it exists
        try:
            exch_file = data_path / '01_EXCH' / 'Network_exchanges.csv'
            if exch_file.exists():
                self.data_reg['Exch']['Network_exchanges'] = pd.read_csv(exch_file, sep=CSV_SEPARATOR,
                                                      header=[0], index_col=[0, 1, 2, 3]).loc[
                                          (self.regions_names, self.regions_names,
                                           self.data_reg['Exch']['Misc_exch']['add_sets']['EXCHANGE_NETWORK_R'],
                                           slice(None)), :]
            else:
                logging.warning('Network_exchanges.csv not found')
        except Exception as e:
            logging.warning(f'Could not read Network_exchanges: {e}')
        
        # Try to read distance data if it exists
        try:
            dist_file = data_path / '01_EXCH' / 'Dist.csv'
            if dist_file.exists():
                self.data_reg['Exch']['Dist'] = pd.read_csv(dist_file, sep=CSV_SEPARATOR,
                                                            header=[0], index_col=[0,1]).loc[
                                                (self.regions_names, self.regions_names), :]
            else:
                logging.warning('Dist.csv not found')
        except Exception as e:
            logging.warning(f'Could not read Dist: {e}')
        # Try to read exchange losses data if it exists
        try:
            exch_losses_file = data_path / '01_EXCH' / 'Exchange_losses.csv'
            if exch_losses_file.exists():
                self.data_reg['Exch']['Exchange_losses'] = pd.read_csv(exch_losses_file, sep=CSV_SEPARATOR,
                                                    header=[0], index_col=[0])
            else:
                logging.warning('Exchange_losses.csv not found')
        except Exception as e:
            logging.warning(f'Could not read Exchange_losses: {e}')
        # Try to read Lhv data if it exists
        try:
            lhv_file = data_path / '01_EXCH' / 'Lhv.csv'
            if lhv_file.exists():
                self.data_reg['Exch']['Lhv'] = pd.read_csv(lhv_file, sep=CSV_SEPARATOR, header=[0], index_col=[0])\
                                           .loc[self.data_reg['Exch']['Misc_exch']['add_sets']['EXCHANGE_FREIGHT_R'], :]
            else:
                logging.warning('Lhv.csv not found')
        except Exception as e:
            logging.warning(f'Could not read Lhv: {e}')        
        return
    
    def concat_reg_data(self, to_concat: list):
        """
        Concatenates across regions the input data
        
        Parameters
        ----------
        to_concat : list
            List of the names of the dataframes to concatenate
        
        Returns
        -------
        tuple
            Tuple containing the concatenated data
        """
        # Create frames for concatenation
        frames = dict.fromkeys(to_concat)
        
        # Misc needs special procedure
        if 'Misc' in to_concat:
            frames['Misc'] = dict.fromkeys(['misc', 'share_ned'])
        
        for c in to_concat:
            if c == 'Misc':
                frames[c]['misc'] = list()
                frames[c]['share_ned'] = list()
            else:
                frames[c] = list()
            
            for n, r in self.regions.items():
                if c == 'Misc':
                    d = r.data[c].copy()
                    share_ned = d.pop('share_ned')
                    frames[c]['misc'].append(pd.Series(d))
                    frames[c]['share_ned'].append(pd.Series(share_ned))
                else:
                    frames[c].append(r.data[c].copy())
        
        # Concatenate and store into a tuple
        out = tuple()
        for c in to_concat:
            if c == 'Misc':
                misc_dict = dict.fromkeys(['misc', 'share_ned'])
                misc_dict['misc'] = pd.concat(
                    frames[c]['misc'],
                    axis=1,
                    keys=self.regions_names,
                    join='inner'
                ).T
                misc_dict['share_ned'] = pd.concat(
                    frames[c]['share_ned'],
                    axis=1,
                    keys=self.regions_names,
                    join='inner'
                ).T
                out = out + (misc_dict,)
            else:
                out = out + (
                    pd.concat(
                        frames[c],
                        axis=0,
                        keys=self.regions_names,
                        join='inner'
                    ),
                )
        
        return out
    
    def build_sets(self):
        """Build all sets from data"""
        logging.info('Building sets from data')
        
        # Demand related sets
        self.sets['SECTORS'] = list(
            self.ref_region.data['Demands']
            .drop(columns=['Category', 'Subcategory', 'Units'], errors='ignore')
            .columns
        )
        self.sets['END_USES_INPUT'] = list(self.ref_region.data['Demands'].index)
        
        # Get END_USES_CATEGORIES from Misc_indep
        if 'add_sets' in self.data_indep.get('Misc_indep', {}):
            if 'END_USES_CATEGORIES' in self.data_indep['Misc_indep']['add_sets']:
                self.sets['END_USES_CATEGORIES'] = list(
                    self.data_indep['Misc_indep']['add_sets']['END_USES_CATEGORIES'].keys()
                )
                self.sets['END_USES_TYPES_OF_CATEGORY'] = \
                    self.data_indep['Misc_indep']['add_sets']['END_USES_CATEGORIES']
            else:
                self.sets['END_USES_CATEGORIES'] = []
                self.sets['END_USES_TYPES_OF_CATEGORY'] = {}
        else:
            self.sets['END_USES_CATEGORIES'] = []
            self.sets['END_USES_TYPES_OF_CATEGORY'] = {}
        
        eut = list(self.sets['END_USES_TYPES_OF_CATEGORY'].values())
        eut = [item for sublist in eut for item in sublist]  # flatten list
        
        # Resources related sets
        self.sets['RESOURCES'] = list(self.data_indep['Resources_indep'].index)
        self.sets['RE_FUELS'] = list(
            self.data_indep['Resources_indep'][
                self.data_indep['Resources_indep'].loc[:, 'Subcategory'] == 'Renewable fuel'
            ].index
        )
        self.sets['RE_RESOURCES'] = list(
            self.data_indep['Resources_indep'][
                self.data_indep['Resources_indep'].loc[:, 'Category'] == 'Renewable'
            ].index
        )
        self.sets['EXPORT'] = list(
            self.data_indep['Resources_indep'][
                self.data_indep['Resources_indep'].loc[:, 'Category'] == 'Export'
            ].index
        )
        
        # Get additional sets from Misc_indep
        if 'add_sets' in self.data_indep.get('Misc_indep', {}):
            add_sets = self.data_indep['Misc_indep']['add_sets']
            self.sets['NOT_LAYERS'] = add_sets.get('NOT_LAYERS', [])
            self.sets['RES_IMPORT_CONSTANT'] = add_sets.get('RES_IMPORT_CONSTANT', [])
            self.sets['EVs_BATT'] = add_sets.get('EVs_BATT', [])
            self.sets['V2G'] = add_sets.get('V2G', [])
            self.sets['STORAGE_DAILY'] = add_sets.get('STORAGE_DAILY', [])
        else:
            self.sets['NOT_LAYERS'] = []
            self.sets['RES_IMPORT_CONSTANT'] = []
            self.sets['EVs_BATT'] = []
            self.sets['V2G'] = []
            self.sets['STORAGE_DAILY'] = []
        
        # Exchanges related sets
        if self.data_reg['Exch'] and 'Misc_exch' in self.data_reg['Exch']:
            exch_add_sets = self.data_reg['Exch']['Misc_exch'].get('add_sets', {})
            self.sets['EXCHANGE_FREIGHT_R'] = exch_add_sets.get('EXCHANGE_FREIGHT_R', [])
            self.sets['EXCHANGE_NETWORK_R'] = exch_add_sets.get('EXCHANGE_NETWORK_R', [])
            self.sets['NETWORK_TYPE'] = exch_add_sets.get('NETWORK_TYPE', [])
        else:
            self.sets['EXCHANGE_FREIGHT_R'] = []
            self.sets['EXCHANGE_NETWORK_R'] = []
            self.sets['NETWORK_TYPE'] = []
        
        self.sets['NOEXCHANGES'] = list(
            set(self.sets['RESOURCES']) -
            set(self.sets['EXCHANGE_FREIGHT_R']) -
            set(self.sets['EXCHANGE_NETWORK_R'])
        )
        
        # Technologies related sets
        all_techs = list(self.ref_region.data['Technologies'].index)
        layers_in_out_tech = self.data_indep['Layers_in_out'].loc[
            ~self.data_indep['Layers_in_out'].index.isin(self.sets['RESOURCES']), :
        ]
        
        self.sets['TECHNOLOGIES_OF_END_USES_TYPE'] = dict.fromkeys(eut)
        for i in eut:
            li = list(layers_in_out_tech.loc[layers_in_out_tech.loc[:, i] == 1, :].index)
            self.sets['TECHNOLOGIES_OF_END_USES_TYPE'][i] = li
        
        all_tech_of_eut = [
            item for sublist in list(self.sets['TECHNOLOGIES_OF_END_USES_TYPE'].values())
            for item in sublist
        ]
        
        self.sets['STORAGE_TECH'] = list(self.data_indep['Storage_eff_in'].index)
        self.sets['INFRASTRUCTURE'] = [
            item for item in all_techs
            if item not in self.sets['STORAGE_TECH']
            and item not in all_tech_of_eut
        ]
        
        self.sets['STORAGE_OF_END_USES_TYPE'] = dict.fromkeys(eut)
        for i in eut:
            li = list(
                self.data_indep['Storage_eff_in'].loc[
                    self.data_indep['Storage_eff_in'].loc[:, i] > 0, :
                ].index
            )
            self.sets['STORAGE_OF_END_USES_TYPE'][i] = li
        
        # Get rid of empty keys
        self.sets['STORAGE_OF_END_USES_TYPE'] = {
            k: v for k, v in self.sets['STORAGE_OF_END_USES_TYPE'].items() if v
        }
        
        # TS_OF_DEC_TECH
        if 'HEAT_LOW_T_DECEN' in self.sets['TECHNOLOGIES_OF_END_USES_TYPE']:
            dec_tech = self.sets['TECHNOLOGIES_OF_END_USES_TYPE']['HEAT_LOW_T_DECEN'].copy()
            if 'DEC_SOLAR' in dec_tech:
                dec_tech.remove('DEC_SOLAR')
            ts_of_dec_tech = [['TS_' + item] for item in dec_tech]
            self.sets['TS_OF_DEC_TECH'] = dict(zip(dec_tech, ts_of_dec_tech))
        else:
            self.sets['TS_OF_DEC_TECH'] = {}
        
        # EVs_BATT_OF_V2G
        self.sets['EVs_BATT_OF_V2G'] = dict(
            zip(self.sets['V2G'], [[item] for item in self.sets['EVs_BATT']])
        )
        
        # BOILERS and COGEN
        self.sets['BOILERS'] = []
        self.sets['COGEN'] = []
        for i in all_tech_of_eut:
            if 'BOILER' in i:
                self.sets['BOILERS'].append(i)
            if 'COGEN' in i:
                self.sets['COGEN'].append(i)
        
        return
    
    def print_indep_dat(self):
        """Print indep.dat file"""
        indep_file = self.output_dir / 'indep.dat'
        logging.info(f'Printing {indep_file}')
        
        # Header
        dp.print_header(
            dat_file=indep_file,
            header_txt=f'File containing data independent of the modelled regions (Year {self.year})'
        )
        
        # Sets
        #dp.newline(
        #    out_path=indep_file,
        #    comment=[
        #        '#----------------------------------------',
        #        '# SETS not depending on TD, nor on REGIONS',
        #        '#----------------------------------------'
        #    ]
        #)
        
        #for n, s in self.sets.items():
        #    if type(s) is list:
        #        dp.print_set(my_set=s, out_path=indep_file, name=n)
        #    else:
        #        for n2, s2 in s.items():
        #            dp.print_set(my_set=s2, out_path=indep_file, name=(n + '["' + n2 + '"]'))
        
        # Parameters
        dp.newline(
            out_path=indep_file,
            comment=[
                '',
                '#----------------------------------------',
                '# PARAMETERS NOT DEPENDING ON THE NUMBER OF TYPICAL DAYS, NOR ON THE REGIONS :',
                '#----------------------------------------'
            ]
        )
        
        for n, d in self.data_indep.items():
            if n != 'Misc_indep':
                df = d.drop(
                    columns=['Category', 'Subcategory', 'Technologies name', 'Units', 'Comment'],
                    errors='ignore'
                )
                df = df.mask(df > 1e14, 'Infinity')
                
                if n == 'Layers_in_out':
                    name = 'param layers_in_out : '
                elif n == 'Storage_eff_in':
                    name = 'param storage_eff_in : '
                elif n == 'Storage_eff_out':
                    name = 'param storage_eff_out : '
                else:
                    name = 'param : '
                
                dp.print_df(
                    df=dp.ampl_syntax(df),
                    out_path=indep_file,
                    name=name,
                    mode='a'
                )
                dp.newline(out_path=indep_file)
            else:
                # Handle Misc_indep
                for n2, d2 in d.items():
                    if type(d2) is dict:
                        if n2 != 'add_sets' and n2 != 'time_series_mapping':
                            if n2 == 'state_of_charge_ev':
                                df = pd.DataFrame.from_dict(
                                    d2,
                                    orient='index',
                                    columns=np.arange(1, 25)
                                )
                                name = 'param state_of_charge_ev : '
                            else:
                                df = pd.DataFrame.from_dict(
                                    d2,
                                    orient='index',
                                    columns=[n2]
                                )
                                name = 'param : '
                            
                            df = df.mask(df > 1e14, 'Infinity')
                            dp.print_df(
                                df=dp.ampl_syntax(df),
                                out_path=indep_file,
                                name=name,
                                mode='a'
                            )
                            dp.newline(out_path=indep_file)
                        else:
                            pass
                    else:
                        if isinstance(d2, (int, float)) and d2 > 1e14:
                            d2 = 'Infinity'
                        dp.print_param(param=d2, out_path=indep_file, name=n2)
                        dp.newline(out_path=indep_file)
        return
    
    def print_reg_data(self):
        """Print regional data files"""
        logging.info('Printing regional data')
        
        # Concatenate data across regions
        self.data_reg['Demands'], \
        self.data_reg['Resources'], \
        self.data_reg['Technologies'], \
        self.data_reg['Storage_power_to_energy'], \
        self.data_reg['Misc'] = self.concat_reg_data(
            to_concat=['Demands', 'Resources', 'Technologies', 'Storage_power_to_energy', 'Misc']
        )
        
        # Print demands, resources, technologies and storage power to energy
        for n in ['Demands', 'Resources', 'Technologies', 'Storage_power_to_energy']:
            self.print_reg_file(n)
        
        # Print misc and exchange files
        self.print_reg_misc()
        self.print_reg_exch()
        
        return
    
    def print_reg_file(self, data_name):
        """Print a regional data file (demands, resources, technologies, storage)"""
        output_file = self.output_dir / f'reg_{data_name.lower()}.dat'
        logging.info(f'Printing {output_file}')
        
        df = self.data_reg[data_name].drop(
            columns=['Category', 'Subcategory', 'Technologies name', 'Units', 'Comment'],
            errors='ignore'
        )
        df = df.mask(df > 1e14, 'Infinity')
        
        if data_name == 'Demands':
            name = 'param end_uses_demand_year : '
        else:
            name = 'param : '
        
        dp.print_df(
            df=dp.ampl_syntax(df),
            out_path=output_file,
            name=name,
            mode='w'
        )
        
        return
    
    def print_reg_misc(self):
        """Print reg_misc.dat file"""
        misc_file = self.output_dir / 'reg_misc.dat'
        logging.info(f'Printing {misc_file}')
        
        # Header
        dp.print_header(
            dat_file=misc_file,
            header_txt='File containing miscellaneous sets and parameters'
        )
        
        # Print REGIONS set
        dp.print_set(my_set=self.regions_names, out_path=misc_file, name='REGIONS')
        dp.newline(out_path=misc_file)
        
        # Get set of regions without dam
        df = self.data_reg['Technologies'].loc[(slice(None), 'HYDRO_DAM'), 'f_max']
        rwithoutdam = list(df.loc[df < 1e-3].index.get_level_values(level=0))
        
        dp.print_set(
            my_set=rwithoutdam,
            out_path=misc_file,
            name='RWITHOUTDAM',
            comment='# Regions without hydro dam'
        )
        dp.newline(out_path=misc_file)
        
        # Print share_ned
        dp.print_df(
            dp.ampl_syntax(self.data_reg['Misc']['share_ned']),
            out_path=misc_file,
            name='param share_ned :'
        )
        
        # Print misc parameters in chunks
        step = 4
        for i in np.arange(0, self.data_reg['Misc']['misc'].shape[1], step):
            df = self.data_reg['Misc']['misc'].iloc[:, i:i + step]
            df = df.mask(df > 1e14, 'Infinity')
            dp.print_df(
                df=dp.ampl_syntax(df),
                out_path=misc_file,
                name='param :'
            )
        
        return
    
    def print_reg_exch(self):
        """Print reg_exch.dat file"""
        exch_file = self.output_dir / 'reg_exch.dat'
        logging.info(f'Printing {exch_file}')
        
        # Header
        dp.print_header(
            dat_file=exch_file,
            header_txt='File containing data related to exchanges between regions'
        )
        
        # Print exchange data
        for n, d in self.data_reg['Exch'].items():
            if n != 'Misc_exch':
                if isinstance(d, pd.DataFrame):
                    dp.print_df(
                        dp.ampl_syntax(d),
                        out_path=exch_file,
                        name='param : '
                    )
            else:
                # Handle Misc_exch
                for n2, d2 in d.items():
                    if type(d2) is dict:
                        if n2 != 'add_sets':
                            df = pd.DataFrame.from_dict(d2, orient='index', columns=[n2])
                            name = 'param : '
                            df = df.mask(df > 1e14, 'Infinity')
                            dp.print_df(
                                df=dp.ampl_syntax(df),
                                out_path=exch_file,
                                name=name,
                                mode='a'
                            )
                            dp.newline(out_path=exch_file)
                    else:
                        if isinstance(d2, (int, float)) and d2 > 1e14:
                            d2 = 'Infinity'
                        dp.print_param(param=d2, out_path=exch_file, name=n2)
                        dp.newline(out_path=exch_file)
        
        return
    
    def generate_all_dat_files(self):
        """Main method to generate all .dat files for this year"""
        logging.info(f'\n{"="*70}')
        logging.info(f'GENERATING .DAT FILES FOR YEAR {self.year}')
        logging.info(f'{"="*70}\n')
        
        try:
            # Read independent data
            self.read_data_indep()
            
            # Initialize regions
            self.init_regions()
            
            # Build sets
            self.build_sets()
            
            # Print indep.dat
            self.print_indep_dat()
            
            # Print regional data files
            self.print_reg_data()
            
            logging.info(f'\n[OK] Successfully generated all .dat files for year {self.year}')
            logging.info(f'Output directory: {self.output_dir}\n')
            
            return True
            
        except Exception as e:
            logging.error(f'\n[ERROR] Failed to generate .dat files for year {self.year}')
            logging.error(f'Error: {e}')
            import traceback
            traceback.print_exc()
            return False


def generate_all_year_dat_files(data_dir=None, regions_names=None, years=None, verbose=True):
    """
    Wrapper function to generate .dat files for all years
    Can be called from other scripts or run standalone
    
    Parameters
    ----------
    data_dir : Path, optional
        Base data directory. If None, uses default path.
    regions_names : list, optional
        List of region names. If None, uses bo_country_code.
    years : list, optional
        List of years to process. If None, uses default years.
    verbose : bool, optional
        Whether to print detailed progress (default: True)
    
    Returns
    -------
    bool
        True if all years processed successfully, False otherwise
    """
    # Set defaults
    if data_dir is None:
        data_dir = Path('/home/pjimenez/ESMC_PTH/EnergyScope_Unif/Data')
    else:
        data_dir = Path(data_dir)
    
    if regions_names is None:
        regions_names = bo_country_code.copy()
    
    if years is None:
        years = ['2015', '2021', '2025', '2030', '2035', '2040', '2045', '2050']
    
    if verbose:
        print(f'\n{"="*70}')
        print('GENERATE INDIVIDUAL YEAR .DAT FILES')
        print(f'{"="*70}\n')
        print(f'Data directory: {data_dir}')
        print(f'Regions: {", ".join(regions_names)}')
        print(f'Years to process: {", ".join(years)}\n')
    else:
        logging.info('Generating individual year .dat files')
        logging.info(f'Data directory: {data_dir}')
        logging.info(f'Years: {", ".join(years)}')
    
    # Check data directory exists
    if not data_dir.exists():
        error_msg = f'Data directory not found: {data_dir}'
        if verbose:
            print(f'[ERROR] {error_msg}')
        else:
            logging.error(error_msg)
        return False
    
    # Process each year
    success_count = 0
    fail_count = 0
    
    for year in years:
        year_dir = data_dir / year
        
        if not year_dir.exists():
            msg = f'Year directory not found: {year_dir}'
            if verbose:
                print(f'[SKIP] {msg}')
            else:
                logging.warning(msg)
            fail_count += 1
            continue
        
        # Create generator and process
        generator = YearDataGenerator(
            year=year,
            regions_names=regions_names,
            data_dir=data_dir
        )
        
        success = generator.generate_all_dat_files()
        
        if success:
            success_count += 1
            
            if verbose:
                # List generated files
                dat_files = list(generator.output_dir.glob('*.dat'))
                print(f'Generated {len(dat_files)} files:')
                for f in sorted(dat_files):
                    print(f'  ✓ {f.name}')
                print()
        else:
            fail_count += 1
    
    # Summary
    if verbose:
        print(f'\n{"="*70}')
        print('SUMMARY')
        print(f'{"="*70}')
        print(f'Successfully processed: {success_count}/{len(years)} years')
        if fail_count > 0:
            print(f'Failed: {fail_count}/{len(years)} years')
        print(f'{"="*70}\n')
        
        if success_count == len(years):
            print('[SUCCESS] All year .dat files generated successfully!')
        elif success_count > 0:
            print('[PARTIAL] Some year .dat files generated successfully')
        else:
            print('[FAILURE] No .dat files were generated')
        
        print('\nFiles generated for each year:')
        print('  1. indep.dat')
        print('  2. reg_demands.dat')
        print('  3. reg_resources.dat')
        print('  4. reg_technologies.dat')
        print('  5. reg_exch.dat')
        print('  6. reg_misc.dat')
        print('  7. reg_storage_power_to_energy.dat')
        print()
    else:
        logging.info(f'Year .dat file generation complete: {success_count}/{len(years)} successful')
        if fail_count > 0:
            logging.warning(f'{fail_count} years failed or skipped')
    
    return success_count == len(years)


def main():
    """Main execution function when run as standalone script"""
    return generate_all_year_dat_files(verbose=True)


if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1)