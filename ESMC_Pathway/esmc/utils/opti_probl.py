"""OptiProbl: set up, solve and read an AMPL optimization problem.

ESMC Pathway - Pablo Jimenez Zabalaga
Based on EnergyScope Multi-Cell (original author: Paolo Thiran).
"""
import numpy as np
import os
import logging
import pandas as pd
import csv
from pathlib import Path
import pickle
from amplpy import AMPL, Environment, DataFrame
from esmc.postprocessing import amplpy2pd as a2p

# TODO: allow choosing the solver together with its options
class OptiProbl:
    """

    The OptiProbl class sets up an optimization problem in AMPL, solves it,
    and interfaces with it through the amplpy API and some additional functions

    Parameters
    ----------
    mod_path : pathlib.Path
        Specifies the path of the .mod file defining the LP problem in ampl syntax
    data_path : list(pathlib.Path)
        List specifying the path of the different .dat files with the data of the LP problem
        in ampl syntax
    options : dict
        Dictionary of the different options for ampl and the cplex solver

    """

    def __init__(self, mod_path=list(), data_path=list(), options=dict(), solver='cplex', ampl_path=None, set_ampl=True):
        """

        Parameters
        ----------
        mod_path
        data_path
        options
        solver
        ampl_path
        set_ampl
        """
        # instantiate different attributes
        if len(mod_path) == 0:
            self.dir = Path()
        else:
            self.dir = mod_path[0].parent

        self.mod_path = mod_path
        self.data_path = data_path
        self.options = options
        if set_ampl:
            self.ampl = self.set_ampl(mod_path, data_path, solver, ampl_path)
        else:
            self.ampl = None
        self.vars = list()
        self.params = list()
        self.sets = dict()
        self.inputs = dict()
        self.t = None
        self.outputs = dict()

        return

    def run_ampl(self):
        """
        Solve the model with AMPL. Writes an IIS file if the model is infeasible.
        """
        try:
            # Solver options
            for o in self.options:
                self.ampl.setOption(o, self.options[o])

            self.ampl.solve()

            # Check solver status
            solve_result = self.ampl.getValue('solve_result_num')
            if solve_result not in [0, 1]:  # 0 = optimal, 1 = feasible
                print("The model is primal/dual infeasible or has no feasible solution.")

                # Try to write the IIS file
                try:
                    print("Trying to write the IIS file...")
                    self.ampl.setOption('cplex_options', 'iisfind=1')
                    self.ampl.eval('write iis;')
                    print("IIS file written to the current working directory.")
                except Exception as e:
                    print(f"Error writing the IIS file: {e}")

            else:
                print("The model was solved successfully.")

            # Reinitialize log options to disable further logging
            self.ampl.setOption('show_stats', 0)
            self.ampl.setOption('times', 0)
            self.ampl.setOption('gentimes', 0)

        except Exception as e:
            print(f"An error occurred during optimization: {e}")
            raise RuntimeError("AMPL optimization failed. See details above.")

    def get_solve_info(self):
        """

       Get the solving info (time and result) and stores it into t attribute

        """
        logging.info('Getting solve_info')
        self.t = list()
        self.t.append(self.ampl.getData('_ampl_elapsed_time;').toList()[0])
        self.t.append(self.ampl.getData('_solve_elapsed_time;').toList()[0])
        self.t.append(self.ampl.getData('solve_result_num;').toList()[0])
        print('[_ampl_elapsed_time, _solve_elapsed_time, solve_result_num]')
        print(self.t)
        # TODO understand why doesn't work with kmedoid_clustering
        return

    def get_inputs(self):
        """

        Get the name of variables and parameters and the sets

        """
        # get values of attributes
        self.get_vars()
        self.get_params()
        self.get_sets()

    def get_vars(self):
        """

        Get the name of the LP optimization problem's variables

        """
        self.vars = list()
        for name, values in self.ampl.getVariables():
            self.vars.append(name)

    def get_params(self):
        """

        Get the name of the LP optimization problem's parameters

               """
        self.params = list()
        for n, p in self.ampl.getParameters():
            self.params.append(n)

    def get_sets(self):
        #TODO update to a more robust version
        """

               Function to sets of the LP optimization problem

        """
        self.sets = dict()
        for name, obj in self.ampl.getSets():
            if len(obj.instances()) <= 1:
                try:
                    self.sets[name] = obj.getValues().toList()
                except Exception as e:
                    logging.warning(str(name) + ' set not working, replacing it by a empty list')
                    self.sets [name] = list()
            else:
                self.sets[name] = self.get_subset(obj)

    def print_inputs(self, directory=None):
        """

        Prints the sets, parameters' names and variables' names of the LP optimization problem

        Parameters
        ----------
        directory : pathlib.Path
        Path of the directory where to save the inputs

        """
        # default directory
        if directory is None:
            directory = self.dir / 'inputs'
        # creating inputs dir
        directory.mkdir(parents=True, exist_ok=True)

        # if params is empty get all inputs
        if not self.params:
            self.get_inputs()
        # printing inputs
        a2p.print_json(self.sets, directory / 'sets.json')
        a2p.print_json(self.params, directory / 'parameters.json')
        a2p.print_json(self.vars, directory / 'variables.json')

        return

    def get_param(self, param_name: str):
        """Function to extract the mentioned parameter and store it into self.inputs

        Parameters
        ----------
        param_name: str
        Name of the parameter to extract from the optimization problem results. Should be written as in the .mod file

        Returns
        -------
        param: pd.DataFrame()
        DataFrame containing the values of the different elements of the parameter.
        The n first columns give the n sets on which it is indexed
        and the last column give the value obtained from the optimization.

        """
        ampl_param = self.ampl.getParameter(param_name)
        # Getting the names of the sets
        indexing_sets = [s.capitalize() for s in ampl_param.getIndexingSets()]
        # Getting the data of the variable into a pandas dataframe
        amplpy_df = ampl_param.getValues()
        param = amplpy_df.toPandas()
        # getting the number of indices. If var has more then 1 index, we set it as a MultiIndex
        # NB: len(indexing_sets) is robust across amplpy versions (getNumIndices() was
        #     removed in newer amplpy releases). Same value, version-independent.
        n_indices = len(indexing_sets)
        if n_indices > 1:
            param.index = pd.MultiIndex.from_tuples(param.index, names=indexing_sets)
        elif n_indices == 1:
            param.index = pd.Index(param.index, name=indexing_sets[0])
        # self.to_pd(ampl_var.getValues()).rename(columns={(var_name+'.val'):var_name})
        self.inputs[param_name] = param
        return param


    def get_var(self, var_name:str):
        """Function to extract the mentioned variable and store it into self.outputs

        Parameters
        ----------
        var_name: str
        Name of the variable to extract from the optimization problem results. Should be written as in the .mod file

        Returns
        -------
        var: pd.DataFrame()
        DataFrame containing the values of the different elements of the variable.
        The n first columns give the n sets on which it is indexed
        and the last column give the value obtained from the optimization.

        """
        ampl_var = self.ampl.getVariable(var_name)
        # Getting the names of the sets
        indexing_sets = [s.capitalize() for s in ampl_var.getIndexingSets()]
        # Getting the data of the variable into a pandas dataframe
        df = ampl_var.get_values().to_pandas()
        df.index.names = indexing_sets
        # getting rid of '.val' (4 trailing characters of the string) into columns names such that the name of the columns correspond to the variable
        df.rename(columns=lambda x: x[:-4], inplace=True)
        #self.to_pd(ampl_var.getValues()).rename(columns={(var_name+'.val'):var_name})
        self.outputs[var_name] = df
        return df

    def read_outputs(self, directory=None):
        """

        Reads the outputs previously printed into csv files to recover a case study without running it again

        Parameters
        ----------
        directory : pathlib.Path
        Path of the directory where the outputs are saved

        """
        # default directory
        if directory is None:
            directory = self.dir / 'outputs'

        with open(directory/'outputs.p', 'rb') as handle:
            self.outputs = pickle.load(handle)
    #############################
    #       STATIC METHODS      #
    #############################

    @staticmethod
    def set_ampl(mod_path=list(), data_path=list(), solver='cplex', ampl_path=None):
        """

        Initialize the AMPL() object containing the LP problem

        Parameters
        ----------
         mod_path : list(pathlib.Path)
        Specifies the path of the .mod files defining the LP problem in ampl syntax

        data_path : list(pathlib.Path)
        List specifying the path of the different .dat files with the data of the LP problem
        in ampl syntax

        solver : str
        Name of the solver (default='cplex')

        ampl_path : None or pathlib.Path
        Default None means ampl is a path variable
        otherwise, give the path to ampl binaries files and solver files

        options : dict
        Dictionary of the different options for ampl and the cplex solver

        Returns
        -------
        ampl object created

        """
        try:
            if ampl_path is None:
                # Create an AMPL instance
                ampl = AMPL()
                # define solver
                ampl.setOption('solver', solver)
            else:
                # Create an AMPL instance
                ampl = AMPL(Environment(binary_directory=str(ampl_path)))#, binary_name='ampl.exe'))
                # define solver
                ampl.setOption('solver', str(ampl_path / solver))

            # Read the model and data files.
            for m in mod_path:
                ampl.read(m)
            for d in data_path:
                ampl.readData(d)
        except Exception as e:
            print(e)
            raise

        return ampl

    @staticmethod
    def get_subset(my_set):
        """

        Function to extract the subsets of set containing sets from the AMPL() object

               Parameters
               ----------
            my_set : amplpy.set.Set
            2-dimensional set to extract


               Returns
               -------
               d : dict()
               dictionary containing the subsets as lists

               """
        d = dict()
        for n, o in my_set.instances():
            try:
                d[n] = o.getValues().toList()
            except Exception as e:
                logging.warning(str(n) + ' subset not working, replacing it by an empty list')
                d[n] = list()
        return d
