#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rolling Horizon Preprocessor for ESMC Pathway Model
====================================================

This module handles the management of optimization windows for rolling horizon
(myopic) pathway optimization in the ESMC model.

Key responsibilities:
- Splitting the full pathway into optimization windows
- Generating PES_seq_opti.dat for each window
- Updating remaining lifetime calculations
- Managing the transition between windows

Based on the AmplPreProcessor from EnergyScope Pathway (single-region),
adapted for the multi-regional ESMC model.

Author: Pablo Jimenez Zabalaga (adapted for ESMC)
Date: January 2026
"""

from pathlib import Path
import os
import re
from copy import deepcopy
import logging
from typing import List, Tuple, Optional


class EsmcPreProcessor:
    """
    Preprocessor for rolling horizon optimization in ESMC pathway model.

    This class manages the sequential optimization of a pathway model by
    splitting it into overlapping time windows. Each window is optimized
    independently, with decisions from previous windows fixed as parameters.

    Parameters
    ----------
    esmc_pathway : EsmcPathway
        The pathway model object containing configuration
    n_years_wnd : int, default=10
        Duration in years of each optimization window (must be multiple of 5)
        Example: 10 means each window covers 10 years (2 phases of 5 years each)
    n_years_overlap : int, default=5
        Duration in years of overlap between consecutive windows (must be multiple of 5)
        Example: 5 means 1 phase overlap between windows
    t_phase : int, default=5
        Duration of each phase in years (typically 5 for pathway models)

    Attributes
    ----------
    years_opti : list of lists
        Years to optimize in each window
    phases_opti : list of lists
        Phases to optimize in each window
    years_up_to : list of lists
        Cumulative years up to and including each window
    phases_up_to : list of lists
        Cumulative phases up to and including each window
    year_to_rm : str
        Year to exclude from results collection (overlap year)

    Example
    -------
    >>> from esmc.utils.esmcpt import EsmcPathway
    >>> from esmc.utils.esmc_preprocessor import EsmcPreProcessor
    >>>
    >>> # Create pathway model
    >>> pathway = EsmcPathway(config, years_config)
    >>>
    >>> # Initialize preprocessor with 10-year windows, 5-year overlap
    >>> pre = EsmcPreProcessor(pathway, n_years_wnd=10, n_years_overlap=5)
    >>>
    >>> # Iterate through windows
    >>> for i in range(pre.get_num_windows()):
    ...     years = pre.write_seq_opti(i)
    ...     pre.update_remaining_years(i)
    ...     # ... optimize and collect results ...

    Notes
    -----
    The window configuration must satisfy:
    - n_years_wnd must be a multiple of t_phase
    - n_years_overlap must be a multiple of t_phase
    - n_years_wnd > n_years_overlap
    - n_years_wnd <= 30 years
    """

    def __init__(self, esmc_pathway, n_years_wnd: int = 10,
                 n_years_overlap: int = 5, t_phase: int = 5):

        # Store reference to model directory
        self.mod_path = esmc_pathway.cs_dir

        # Extract years from pathway model
        # Format: ['YEAR_2015', 'YEAR_2021', 'YEAR_2025', ...]
        self.list_year = ['YEAR_' + y for y in esmc_pathway.years]

        # Extract phases from pathway model
        # Format: ['2015_2021', '2021_2025', '2025_2030', ...]
        self.list_phase = list(esmc_pathway.phases.keys())

        # Store original lists before modification
        self.list_year_full = self.list_year.copy()
        self.list_phase_full = self.list_phase.copy()

        # Remove initial year/phase (2015 is fixed initial condition)
        self._remove_initial_year_phase()

        # Store window configuration
        self.t_phase = t_phase
        self.n_years_wnd = n_years_wnd
        self.n_years_overlap = n_years_overlap

        # Initialize window definition lists
        self.years_opti: List[List[str]] = []
        self.phases_opti: List[List[str]] = []
        self.years_up_to: List[List[str]] = []
        self.phases_up_to: List[List[str]] = []
        self.year_to_rm: str = ''

        # Store regions reference
        self.regions = esmc_pathway.regions_names

        # Validate configuration and compute windows
        self._validate_window_config()
        self._compute_pathway_windows()

        # Log initialization
        logging.info(f'EsmcPreProcessor initialized:')
        logging.info(f'  Total years (excluding 2015): {len(self.list_year)}')
        logging.info(f'  Total phases (excluding 2015_2021): {len(self.list_phase)}')
        logging.info(f'  Window size: {n_years_wnd} years')
        logging.info(f'  Overlap: {n_years_overlap} years')
        logging.info(f'  Number of windows: {len(self.years_opti)}')

    def _remove_initial_year_phase(self) -> None:
        """
        Remove YEAR_2015 and 2015_2021 from optimization lists.

        The year 2015 and phase 2015_2021 represent fixed initial conditions
        and should not be part of the optimization.
        """
        if 'YEAR_2015' in self.list_year:
            self.list_year.remove('YEAR_2015')
            logging.debug('Removed YEAR_2015 from optimization years')
        if '2015_2021' in self.list_phase:
            self.list_phase.remove('2015_2021')
            logging.debug('Removed 2015_2021 from optimization phases')

    def _validate_window_config(self) -> None:
        """
        Validate window configuration parameters.

        Raises
        ------
        ValueError
            If any configuration parameter is invalid
        """
        errors = []

        if self.n_years_wnd % self.t_phase != 0:
            errors.append(f'Window size ({self.n_years_wnd}) must be multiple of t_phase ({self.t_phase})')

        if self.n_years_overlap % self.t_phase != 0:
            errors.append(f'Overlap ({self.n_years_overlap}) must be multiple of t_phase ({self.t_phase})')

        if self.n_years_wnd > 30 or self.n_years_wnd <= 0:
            errors.append(f'Window size must be between 1 and 30 years (got {self.n_years_wnd})')

        if self.n_years_wnd < self.t_phase:
            errors.append(f'Window size ({self.n_years_wnd}) must be >= t_phase ({self.t_phase})')

        if self.n_years_wnd <= self.n_years_overlap:
            errors.append(f'Window size ({self.n_years_wnd}) must be > overlap ({self.n_years_overlap})')

        if errors:
            raise ValueError('Invalid window configuration:\n  ' + '\n  '.join(errors))

        logging.debug('Window configuration validated successfully')

    def _compute_pathway_windows(self) -> None:
        """
        Split the pathway into optimization windows.

        This method computes:
        - years_opti: Which years belong to each window
        - phases_opti: Which phases belong to each window
        - years_up_to: Cumulative years up to each window (for fixing)
        - phases_up_to: Cumulative phases up to each window (for fixing)
        """
        # Calculate number of phases/years per window and overlap
        n_phase_window = self.n_years_wnd // self.t_phase
        n_year_window = self.n_years_wnd // self.t_phase
        n_phase_overlap = self.n_years_overlap // self.t_phase
        n_year_overlap = self.n_years_overlap // self.t_phase

        logging.debug(f'Phases per window: {n_phase_window}, overlap: {n_phase_overlap}')

        # Initialize lists with maximum possible size
        max_windows = len(self.list_year)
        years_opti = [[] for _ in range(max_windows)]
        phases_opti = [[] for _ in range(max_windows)]

        # Compute years for each window
        for i in range(max_windows):
            start_idx = i * n_year_window - i * n_year_overlap
            end_idx = (n_year_window + 1) * (i + 1) - i * (n_year_overlap + 1)

            if start_idx >= len(self.list_year):
                break

            years_opti[i] = self.list_year[start_idx:min(end_idx, len(self.list_year))]

            if end_idx >= len(self.list_year):
                break

        # Compute phases for each window
        for i in range(max_windows):
            start_idx = i * n_phase_window - i * n_phase_overlap
            end_idx = n_phase_window * (i + 1) - i * n_phase_overlap

            if start_idx >= len(self.list_phase):
                break

            phases_opti[i] = self.list_phase[start_idx:min(end_idx, len(self.list_phase))]

            if end_idx >= len(self.list_phase):
                break

        # Remove empty windows
        years_opti = [ele for ele in years_opti if ele]
        phases_opti = [ele for ele in phases_opti if ele]

        # Compute cumulative "up to" lists
        years_up_to = [[] for _ in range(len(years_opti))]
        phases_up_to = [[] for _ in range(len(phases_opti))]

        for i in range(len(years_up_to) - 1):
            end_idx = (n_year_window + 1) * (i + 1) - i * (n_year_overlap + 1) - n_year_overlap
            years_up_to[i] = self.list_year[0:end_idx]
        years_up_to[-1] = self.list_year.copy()

        for i in range(len(phases_up_to) - 1):
            end_idx = n_phase_window * (i + 1) - i * n_phase_overlap - n_phase_overlap
            phases_up_to[i] = self.list_phase[0:end_idx]
        phases_up_to[-1] = self.list_phase.copy()

        # Store results
        self.years_opti = years_opti
        self.phases_opti = phases_opti
        self.years_up_to = years_up_to
        self.phases_up_to = phases_up_to

        # Log window configuration
        logging.info('Window configuration computed:')
        for i, (years, phases) in enumerate(zip(self.years_opti, self.phases_opti)):
            logging.info(f'  Window {i+1}: {len(years)} years ({years[0]} to {years[-1]}), '
                        f'{len(phases)} phases')

    def write_seq_opti(self, window_idx: int) -> List[str]:
        """
        Write PES_seq_opti.dat for the specified optimization window.

        This file defines the AMPL sets and parameters for the optimization:
        - YEARS, PHASE: Full sets for all years/phases
        - PHASE_START, PHASE_STOP: Start/stop years for each phase
        - YEARS_WND, PHASE_WND: Window-specific sets
        - YEARS_UP_TO, PHASE_UP_TO: Cumulative sets
        - YEAR_ONE, YEAR_ONE_NEXT: For fixing variables
        - t_phase, diff_2015_phase, phase_to_year: Parameters

        Parameters
        ----------
        window_idx : int
            Index of current window (0-based)

        Returns
        -------
        list of str
            Years included in the current optimization window
        """
        if window_idx < 0 or window_idx >= len(self.years_opti):
            raise ValueError(f'Invalid window index: {window_idx}. '
                           f'Valid range: 0-{len(self.years_opti)-1}')

        curr_years_wnd = self.years_opti[window_idx]
        curr_phases_wnd = self.phases_opti[window_idx]
        curr_years_up_to = self.years_up_to[window_idx]
        curr_phases_up_to = self.phases_up_to[window_idx]

        year_one = self.years_opti[window_idx][0]

        # Track year to remove for results collection (overlap year)
        if window_idx > 0:
            self.year_to_rm = year_one

        # Determine next window's first year (for overlap)
        if window_idx == len(self.years_opti) - 1:
            has_next_window = False
            year_one_next = ''
        else:
            has_next_window = True
            year_one_next = self.years_opti[window_idx + 1][0]

        # Write the .dat file
        output_file = self.mod_path / 'PES_seq_opti.dat'

        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(f'# Sequential optimization data for window {window_idx + 1}\n')
            f.write(f'# Generated by EsmcPreProcessor\n')
            f.write(f'# Years in window: {curr_years_wnd}\n')
            f.write(f'# Phases in window: {curr_phases_wnd}\n\n')

            # =====================================================
            # FULL SETS (needed for all data files to load correctly)
            # =====================================================

            # Write full YEARS set (including YEAR_2015)
            f.write('# Years in the analysis\n')
            f.write('set YEARS := ')
            f.write(' '.join(self.list_year_full))
            f.write(';\n\n')

            # Write full PHASE set (including 2015_2021)
            f.write('# Optimization phases\n')
            f.write('set PHASE := ')
            f.write(' '.join(self.list_phase_full))
            f.write(';\n\n')

            # Write PHASE_START and PHASE_STOP for each phase
            f.write('# Phase start and stop years\n')
            for phase in self.list_phase_full:
                years_in_phase = phase.split('_')
                start_year = f'YEAR_{years_in_phase[0]}'
                stop_year = f'YEAR_{years_in_phase[1]}'
                f.write(f'set PHASE_START[{phase}] := {start_year};\n')
                f.write(f'set PHASE_STOP[{phase}] := {stop_year};\n\n')

            # =====================================================
            # WINDOW-SPECIFIC SETS
            # =====================================================
            f.write('# Window-specific sets\n')

            # Write YEARS_WND set
            f.write('set YEARS_WND := ')
            f.write(' '.join(curr_years_wnd))
            f.write(';\n\n')

            # Write PHASE_WND set
            f.write('set PHASE_WND := ')
            f.write(' '.join(curr_phases_wnd))
            f.write(';\n\n')

            # Write YEARS_UP_TO set
            f.write('set YEARS_UP_TO := ')
            f.write(' '.join(curr_years_up_to))
            f.write(';\n\n')

            # Write PHASE_UP_TO set
            f.write('set PHASE_UP_TO := ')
            f.write(' '.join(curr_phases_up_to))
            f.write(';\n\n')

            # Write YEAR_ONE set
            # Window 1: YEAR_2015 (fixed initial condition)
            # Windows 2+: First year of current window (to be fixed)
            f.write('set YEAR_ONE')
            if window_idx == 0:
                f.write(' := YEAR_2015')
            else:
                f.write(f' := {year_one}')
            f.write(';\n\n')

            # Write YEAR_ONE_NEXT set (empty for last window)
            f.write('set YEAR_ONE_NEXT')
            if has_next_window:
                f.write(f' := {year_one_next}')
            f.write(';\n\n')

            # =====================================================
            # PARAMETERS
            # =====================================================

            # Write t_phase parameter
            f.write('# Phase durations [years]\n')
            f.write('param t_phase :=\n')
            phase_durations = {
                '2015_2021': 6,
                '2021_2025': 4,
                '2025_2030': 5,
                '2030_2035': 5,
                '2035_2040': 5,
                '2040_2045': 5,
                '2045_2050': 5
            }
            for phase in self.list_phase_full:
                duration = phase_durations.get(phase, 5)
                f.write(f'{phase}    {duration}\n')
            f.write(';\n\n')

            # Write diff_2015_phase parameter
            f.write('# Years from 2015 to each phase (for discounting)\n')
            f.write('param diff_2015_phase :=\n')
            diff_values = {
                '2015_2021': 3,      # midpoint: 2018 - 2015 = 3
                '2021_2025': 8,      # midpoint: 2023 - 2015 = 8
                '2025_2030': 12.5,   # midpoint: 2027.5 - 2015 = 12.5
                '2030_2035': 17.5,   # midpoint: 2032.5 - 2015 = 17.5
                '2035_2040': 22.5,   # midpoint: 2037.5 - 2015 = 22.5
                '2040_2045': 27.5,   # midpoint: 2042.5 - 2015 = 27.5
                '2045_2050': 32.5    # midpoint: 2047.5 - 2015 = 32.5
            }
            for phase in self.list_phase_full:
                diff = diff_values.get(phase, 0)
                f.write(f'{phase}    {diff}\n')
            f.write(';\n\n')

            # Write phase_to_year parameter
            f.write('# Mapping phase to year\n')
            f.write('param phase_to_year :=\n')
            for phase in self.list_phase_full:
                end_year = phase.split('_')[1]
                f.write(f'{phase}    {end_year}\n')
            f.write(';\n')

        logging.info(f'Wrote PES_seq_opti.dat for window {window_idx + 1}:')
        logging.info(f'  YEARS_WND: {curr_years_wnd}')
        logging.info(f'  PHASE_WND: {curr_phases_wnd}')

        return curr_years_wnd.copy()

        return curr_years_wnd.copy()

    def update_remaining_years(self, window_idx: int,
                               file_in: str = 'PES_data_remaining.dat',
                               file_out: str = 'PES_data_remaining_wnd.dat') -> None:
        """
        Update remaining lifetime calculations for current optimization window.

        The remaining years parameter determines how much investment return
        is credited for technologies that survive past the end of the horizon
        (typically 2050).

        Parameters
        ----------
        window_idx : int
            Index of current window (0-based)
        file_in : str, default='PES_data_remaining.dat'
            Input file with full remaining years data for all phases
        file_out : str, default='PES_data_remaining_wnd.dat'
            Output file with window-specific remaining years data
        """
        # Build list of phases to include
        curr_phases_up_to = ['2015_2021'] + self.phases_up_to[window_idx]
        curr_phases_wnd = self.phases_opti[window_idx]
        phase_list = deepcopy(curr_phases_up_to)

        last_phase = self.list_phase[-1]

        # Add overlap phases if not last window
        if self.n_years_overlap > 0 and curr_phases_wnd[-1] != last_phase:
            n_phase_extra = self.n_years_overlap // self.t_phase
            phase_list += curr_phases_wnd[-n_phase_extra:]

        n_phase = len(phase_list)

        # File paths
        pth_file_in = self.mod_path / file_in
        pth_file_out = self.mod_path / file_out

        # Check if input file exists
        if not pth_file_in.exists():
            logging.warning(f'Remaining years file not found: {pth_file_in}')
            logging.warning('Creating empty output file')
            with open(pth_file_out, 'w') as f:
                f.write(f'# No remaining years data available\n')
            return

        # Read original file structure
        with open(pth_file_in, encoding='utf-8') as fp:
            header_line = next(fp)
            lines = fp.readlines()
            if lines:
                n_cols = len(re.split(r'\t+', lines[0].rstrip('\t')))
            else:
                n_cols = 0

        # Write window-specific file
        with open(pth_file_out, 'w', encoding='utf-8') as f:
            f.write('param remaining_years : ')
            f.write(' '.join(phase_list))
            f.write(' := \n')

            with open(pth_file_in, encoding='utf-8') as fp:
                next(fp)  # Skip header
                for line in fp:
                    parts = re.split(r'\t+', line.rstrip('\t'))
                    if len(parts) < n_cols:
                        break

                    # Write technology name
                    f.write(f'{parts[0]}\t')

                    # Write relevant phase columns
                    for i in range(len(parts) - n_phase - 1, len(parts) - 1):
                        f.write(f'{parts[i]}\t')
                    f.write('\n')

            f.write(';\n')

        logging.info(f'Updated remaining years for window {window_idx + 1}')
        logging.debug(f'  Phases included: {phase_list}')

    def get_num_windows(self) -> int:
        """
        Get the total number of optimization windows.

        Returns
        -------
        int
            Number of windows in the rolling horizon
        """
        return len(self.years_opti)

    def get_window_info(self, window_idx: int) -> dict:
        """
        Get detailed information about a specific window.

        Parameters
        ----------
        window_idx : int
            Index of window (0-based)

        Returns
        -------
        dict
            Dictionary with window information:
            - 'years': Years in window
            - 'phases': Phases in window
            - 'years_up_to': Cumulative years
            - 'phases_up_to': Cumulative phases
            - 'is_first': Whether this is the first window
            - 'is_last': Whether this is the last window
        """
        if window_idx < 0 or window_idx >= len(self.years_opti):
            raise ValueError(f'Invalid window index: {window_idx}')

        return {
            'index': window_idx,
            'years': self.years_opti[window_idx],
            'phases': self.phases_opti[window_idx],
            'years_up_to': self.years_up_to[window_idx],
            'phases_up_to': self.phases_up_to[window_idx],
            'is_first': window_idx == 0,
            'is_last': window_idx == len(self.years_opti) - 1,
            'n_years': len(self.years_opti[window_idx]),
            'n_phases': len(self.phases_opti[window_idx])
        }

    def print_window_summary(self) -> None:
        """Print a summary of all optimization windows."""
        print(f"\n{'='*60}")
        print(f"Rolling Horizon Window Summary")
        print(f"{'='*60}")
        print(f"Configuration:")
        print(f"  Window size: {self.n_years_wnd} years")
        print(f"  Overlap: {self.n_years_overlap} years")
        print(f"  Phase duration: {self.t_phase} years")
        print(f"\nWindows: {self.get_num_windows()}")
        print(f"{'-'*60}")

        for i in range(self.get_num_windows()):
            info = self.get_window_info(i)
            marker = " [FIRST]" if info['is_first'] else (" [LAST]" if info['is_last'] else "")
            print(f"\nWindow {i+1}{marker}:")
            print(f"  Years: {info['years'][0]} to {info['years'][-1]} ({info['n_years']} years)")
            print(f"  Phases: {info['phases'][0]} to {info['phases'][-1]} ({info['n_phases']} phases)")
            print(f"  Years up to: {len(info['years_up_to'])} years cumulative")

        print(f"\n{'='*60}")