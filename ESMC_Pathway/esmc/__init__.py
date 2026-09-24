"""ESMC Pathway: multi-regional energy system pathway model (EnergyScope Multi-Cell).

Importing this package sets up logging:
- console: INFO messages
- file: one log per run in <repository root>/logs/ (name from commons['logfile'])
"""
import logging.config
from pathlib import Path

from .common import commons, PROJECT_DIR

LOG_DIR = PROJECT_DIR / 'logs'
try:
    LOG_DIR.mkdir(exist_ok=True)
except OSError:
    LOG_DIR = Path.cwd()   # read-only install: write logs in the working directory

logging.config.dictConfig({
    'version': 1,
    'disable_existing_loggers': False,
    'formatters': {
        'standard': {'format': '%(asctime)s [%(levelname)-8s] (%(funcName)s): %(message)s',
                     'datefmt': '%y/%m/%d %H:%M:%S'},
        'notime': {'format': '[%(levelname)-8s] (%(funcName)s): %(message)s'},
    },
    'handlers': {
        'console': {'class': 'logging.StreamHandler', 'stream': 'ext://sys.stderr',
                    'level': 'INFO', 'formatter': 'notime'},
        'file': {'class': 'logging.FileHandler', 'level': 'INFO', 'formatter': 'standard',
                 'filename': str(LOG_DIR / commons['logfile']), 'encoding': 'utf8',
                 'delay': True},   # file created only when something is logged
    },
    'root': {'level': 'INFO', 'handlers': ['console', 'file']},
})
