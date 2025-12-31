"""
AMPL Results Cache System (Updated with Cost Phase Variables)
==============================================================

This module provides functionality to save and load AMPL optimization results,
allowing you to skip the 2.5-hour optimization and directly run post-processing.

Updated to include C_inv_phase_tech, C_op_phase_tech, and C_op_phase_res for cost plotting.

Usage:
------
1. After optimization, call: save_ampl_results(pathway_model)
2. To skip optimization, call: load_ampl_results(pathway_model)

This saves ~2.5 hours when iterating on result processing/printing code.
"""

import pickle
import logging
from pathlib import Path
import pandas as pd
from datetime import datetime


def save_ampl_results(pathway_model, cache_file=None):
    """
    Save all AMPL variables and parameters to a cache file
    
    This function extracts and saves all variables and parameters from AMPL
    that are needed for post-processing, so you can skip re-optimization.
    
    Parameters:
    -----------
    pathway_model : EsmcPathway
        The pathway model after optimization
    cache_file : str or Path, optional
        Path to save the cache. If None, saves to case_study directory
        
    Returns:
    --------
    Path : Location where cache was saved
    """
    logging.info("Saving AMPL results to cache...")
    
    # Determine cache location
    if cache_file is None:
        cache_file = pathway_model.cs_dir / 'ampl_results_cache.pkl'
    else:
        cache_file = Path(cache_file)
    
    # Create backup directory if it doesn't exist
    cache_file.parent.mkdir(parents=True, exist_ok=True)
    
    # Dictionary to store all results
    cache = {
        'metadata': {
            'timestamp': datetime.now().isoformat(),
            'case_study': pathway_model.case_study,
            'nbr_td': getattr(pathway_model, 'nbr_td', getattr(pathway_model, 'Nbr_TD', 16)),
            'regions_names': pathway_model.regions_names if hasattr(pathway_model, 'regions_names') else None,
        },
        'variables': {},
        'parameters': {},
        'solve_info': {}
    }
    
    # Get AMPL object
    ampl = pathway_model.esom.ampl
    
    # Save solver status
    try:
        cache['solve_info']['solve_result'] = ampl.getValue("solve_result")
        logging.info(f"  Solver status: {cache['solve_info']['solve_result']}")
    except:
        cache['solve_info']['solve_result'] = 'unknown'
    
    # List of variables to extract (INCLUDING NEW COST PHASE VARIABLES)
    variables_to_save = [
        'F',                    # Installed capacity
        'F_t',                  # Hourly production
        'TotalCost',           # Total cost
        'C_inv',               # Investment cost
        'C_maint',             # Maintenance cost
        'C_op',                # Operation cost
        'GWP_constr',          # GWP from construction
        'GWP_op',              # GWP from operation
        'CO2_net',             # Net CO2 emissions
        'R_t_local',           # Local resources
        'R_t_exterior',        # Exterior resources
        'Exch_imp',            # Imports exchanges
        'Exch_exp',            # Exports exchanges
        'Transfer_capacity',   # Transfer capacity
        'Storage_in',          # Storage input
        'Storage_out',         # Storage output
        'Storage_level',       # Storage level
        'End_uses',            # End uses
        'Curt',                # Curtailment
        'F_new',               # New capacity (pathway)
        'F_old',               # Retired capacity (pathway)
        'F_ant_old',           # Retired capacity from anterior phases (pathway)
        'F_decom',             # Decommissioned capacity (pathway)
        'C_inv_phase_tech',    # Investment per tech/region/phase (pathway) - NEW
        'C_op_phase_tech',     # Operation cost per tech/region/phase (pathway) - NEW
        'C_op_phase_res',      # Resource cost per region/phase (pathway) - NEW
    ]
    
    # List of parameters to extract
    parameters_to_save = [
        't_op',                # Operating time
        'tau',                 # Annualization factor
        'lifetime',            # Technology lifetime
        'layers_in_out',       # Technology layers
    ]
    
    # Extract variables
    logging.info("  Extracting variables...")
    for var_name in variables_to_save:
        try:
            var = pathway_model.esom.get_var(var_name)
            if var is not None and not var.empty:
                cache['variables'][var_name] = var
                logging.info(f"    ✓ {var_name}: {len(var)} rows")
            else:
                logging.warning(f"    ⚠ {var_name}: empty or None")
        except Exception as e:
            logging.warning(f"    ⚠ {var_name}: not found ({e})")
    
    # Extract parameters
    logging.info("  Extracting parameters...")
    for param_name in parameters_to_save:
        try:
            param = pathway_model.esom.get_param(param_name)
            if param is not None and not param.empty:
                cache['parameters'][param_name] = param
                logging.info(f"    ✓ {param_name}: {len(param)} rows")
            else:
                logging.warning(f"    ⚠ {param_name}: empty or None")
        except Exception as e:
            logging.warning(f"    ⚠ {param_name}: not found ({e})")
    
    # Save td_of_days if available
    if hasattr(pathway_model, 'td_of_days'):
        cache['td_of_days'] = pathway_model.td_of_days
        logging.info(f"    ✓ td_of_days saved")
    
    # Save t_h_td if available
    if hasattr(pathway_model, 't_h_td'):
        cache['t_h_td'] = pathway_model.t_h_td
        logging.info(f"    ✓ t_h_td saved")
    
    # Save to pickle file
    try:
        with open(cache_file, 'wb') as f:
            pickle.dump(cache, f, protocol=pickle.HIGHEST_PROTOCOL)
        
        file_size_mb = cache_file.stat().st_size / (1024 * 1024)
        logging.info(f"\n✓ AMPL results cached successfully!")
        logging.info(f"  Location: {cache_file}")
        logging.info(f"  Size: {file_size_mb:.1f} MB")
        logging.info(f"  Variables: {len(cache['variables'])}")
        logging.info(f"  Parameters: {len(cache['parameters'])}")
        
        return cache_file
        
    except Exception as e:
        logging.error(f"✗ Failed to save cache: {e}")
        raise


def load_ampl_results(pathway_model, cache_file=None):
    """
    Load AMPL variables and parameters from cache file
    
    This function loads previously saved AMPL results, allowing you to
    skip the 2.5-hour optimization and go directly to post-processing.
    
    Parameters:
    -----------
    pathway_model : EsmcPathway
        The pathway model (must be initialized but not optimized)
    cache_file : str or Path, optional
        Path to the cache file. If None, looks in case_study directory
        
    Returns:
    --------
    bool : True if successfully loaded, False otherwise
    """
    logging.info("Loading AMPL results from cache...")
    
    # Determine cache location
    if cache_file is None:
        cache_file = pathway_model.cs_dir / 'ampl_results_cache.pkl'
    else:
        cache_file = Path(cache_file)
    
    # Check if file exists
    if not cache_file.exists():
        logging.error(f"✗ Cache file not found: {cache_file}")
        logging.error("  Run optimization first with save_ampl_results()")
        return False
    
    # Load cache
    try:
        with open(cache_file, 'rb') as f:
            cache = pickle.load(f)
        
        file_size_mb = cache_file.stat().st_size / (1024 * 1024)
        logging.info(f"  ✓ Loaded cache file: {cache_file}")
        logging.info(f"    Size: {file_size_mb:.1f} MB")
        logging.info(f"    Created: {cache['metadata']['timestamp']}")
        logging.info(f"    Case: {cache['metadata']['case_study']}")
        
    except Exception as e:
        logging.error(f"✗ Failed to load cache: {e}")
        return False
    
    # Create a mock AMPL interface that returns cached data
    class CachedAMPLInterface:
        def __init__(self, cached_variables, cached_parameters, solve_info):
            self.cached_variables = cached_variables
            self.cached_parameters = cached_parameters
            self.solve_info = solve_info
        
        def get_var(self, var_name):
            """Get a cached variable"""
            if var_name in self.cached_variables:
                return self.cached_variables[var_name].copy()
            else:
                logging.warning(f"Variable {var_name} not in cache")
                return pd.DataFrame()
        
        def get_param(self, param_name):
            """Get a cached parameter"""
            if param_name in self.cached_parameters:
                return self.cached_parameters[param_name].copy()
            else:
                logging.warning(f"Parameter {param_name} not in cache")
                return pd.DataFrame()
        
        def getValue(self, name):
            """Get a cached value"""
            if name == "solve_result":
                return self.solve_info.get('solve_result', 'cached')
            return None
    
    # Replace the AMPL interface with cached version
    if not hasattr(pathway_model, 'esom') or pathway_model.esom is None:
        # Create a minimal esom object
        class MinimalEsom:
            pass
        pathway_model.esom = MinimalEsom()
    
    # Inject the cached interface
    pathway_model.esom.get_var = CachedAMPLInterface(
        cache['variables'], 
        cache['parameters'],
        cache['solve_info']
    ).get_var
    
    pathway_model.esom.get_param = CachedAMPLInterface(
        cache['variables'], 
        cache['parameters'],
        cache['solve_info']
    ).get_param
    
    # Create a mock ampl object for getValue
    # Create mock classes that mimic the full AMPL API chain
    class MockAMPLDataFrame:
        """Mock AMPL DataFrame that wraps a pandas DataFrame"""
        def __init__(self, pandas_df):
            self.pandas_df = pandas_df if pandas_df is not None else pd.DataFrame()
        
        def toPandas(self):
            """Convert to pandas DataFrame (already is one, just return it)"""
            return self.pandas_df
        
        def toDict(self):
            """Convert to dictionary"""
            return self.pandas_df.to_dict()
    
    class MockVariable:
        def __init__(self, data, var_name=''):
            self.data = data if data is not None else pd.DataFrame()
            self.var_name = var_name
        
        def getValues(self):
            df = self.data.copy()
            if not df.empty and self.var_name:
                column_mapping = {}
                for col in df.columns:
                    if col not in df.index.names:
                        column_mapping[col] = f"{col}.val"  # Add .val suffix
                if column_mapping:
                    df = df.rename(columns=column_mapping)
            return MockAMPLDataFrame(df)
    
    class MockParameter:
        def __init__(self, data, param_name=''):
            self.data = data if data is not None else pd.DataFrame()
            self.param_name = param_name
        
        def getValues(self):
            df = self.data.copy()
            if not df.empty and self.param_name:
                column_mapping = {}
                for col in df.columns:
                    if col not in df.index.names:
                        column_mapping[col] = f"{col}.val"  # Add .val suffix
                if column_mapping:
                    df = df.rename(columns=column_mapping)
            return MockAMPLDataFrame(df)
    
    class MockAmpl:
        def __init__(self, solve_info, cached_variables, cached_parameters):
            self.solve_info = solve_info
            self.cached_variables = cached_variables
            self.cached_parameters = cached_parameters
        
        def getValue(self, name):
            """Get a cached value"""
            if name == "solve_result":
                return self.solve_info.get('solve_result', 'cached')
            return None
        
        def getVariable(self, var_name):
            if var_name in self.cached_variables:
                return MockVariable(self.cached_variables[var_name].copy(), var_name)
            else:
                return MockVariable(pd.DataFrame(), var_name)
        
        def getParameter(self, param_name):
            """Get a cached parameter (AMPL API compatibility)"""
            if param_name in self.cached_parameters:
                return MockParameter(self.cached_parameters[param_name].copy(), param_name)
            else:
                logging.warning(f"Parameter {param_name} not in cache, returning empty")
                return MockParameter(pd.DataFrame(), param_name)
        
        def getSet(self, set_name):
            """Get a set from AMPL (mock - returns empty for cached mode)"""
            # Sets are not typically cached, so return a mock empty set
            logging.debug(f"getSet({set_name}) called on MockAmpl - returning empty")
            
            class MockSet:
                def getValues(self):
                    return MockAMPLDataFrame(pd.DataFrame())
                def toList(self):
                    return []
            
            return MockSet()
        
        def close(self):
            """Mock close - does nothing"""
            pass
    
    pathway_model.esom.ampl = MockAmpl(cache['solve_info'], cache['variables'], cache['parameters'])
    
    # Restore td_of_days and t_h_td if available
    if 'td_of_days' in cache:
        pathway_model.td_of_days = cache['td_of_days']
        logging.info("  ✓ Restored td_of_days")
    
    if 't_h_td' in cache:
        pathway_model.t_h_td = cache['t_h_td']
        logging.info("  ✓ Restored t_h_td")
    
    logging.info(f"\n✓ Successfully loaded {len(cache['variables'])} variables and {len(cache['parameters'])} parameters")
    logging.info("  You can now run post-processing without re-optimizing")
    
    return True


def check_cache_exists(pathway_model, cache_file=None):
    """
    Check if a cache file exists for this model
    
    Parameters:
    -----------
    pathway_model : EsmcPathway
        The pathway model
    cache_file : str or Path, optional
        Path to check. If None, checks default location
        
    Returns:
    --------
    tuple : (exists: bool, file_path: Path, file_info: dict)
    """
    if cache_file is None:
        cache_file = pathway_model.cs_dir / 'ampl_results_cache.pkl'
    else:
        cache_file = Path(cache_file)
    
    exists = cache_file.exists()
    
    file_info = {}
    if exists:
        file_info['size_mb'] = cache_file.stat().st_size / (1024 * 1024)
        file_info['modified'] = datetime.fromtimestamp(cache_file.stat().st_mtime)
        
        # Try to read metadata
        try:
            with open(cache_file, 'rb') as f:
                cache = pickle.load(f)
            file_info['metadata'] = cache.get('metadata', {})
            file_info['n_variables'] = len(cache.get('variables', {}))
            file_info['n_parameters'] = len(cache.get('parameters', {}))
        except:
            pass
    
    return exists, cache_file, file_info


# Example usage in run_pathway.py:
"""
# After optimization (Step 5):
if optimization_successful:
    cache_file = save_ampl_results(pathway_model)
    print(f"Results cached to: {cache_file}")
    print("Next time, you can skip optimization by using --skip-optimize flag")

# To skip optimization:
# python run_pathway.py --skip-optimize

# In main():
if args.skip_optimize:
    if load_ampl_results(pathway_model):
        print("Loaded cached results, skipping optimization")
        # Jump directly to Step 7 (extraction)
    else:
        print("No cache found, must run optimization")
        return
"""