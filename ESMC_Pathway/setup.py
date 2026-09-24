from pathlib import Path
from setuptools import setup, find_packages

HERE = Path(__file__).resolve().parent

setup(
    name='esmc_pathway',
    version=(HERE / 'VERSION').read_text().strip(),
    author='Pablo Jimenez Zabalaga',
    description='ESMC Pathway: multi-regional energy transition pathway model based on EnergyScope Multi-Cell',
    long_description=(HERE / 'README.md').read_text(encoding='utf-8'),
    long_description_content_type='text/markdown',
    license='Apache-2.0',
    packages=find_packages(),
    py_modules=['run'],
    include_package_data=True,
    package_data={'esmc': ['energy_model/*.mod', 'energy_model/headers/*.txt',
                           'preprocessing/kmedoid_clustering/*']},
    python_requires='>=3.9',
    install_requires=[
        'numpy',
        'pandas',
        'amplpy',
        'plotly',
        'kaleido',
    ],
    extras_require={'segmentation': ['tsam']},
)
