from setuptools import setup, find_packages

VERSION = '0.2'

DESCRIPTION = 'BitBirch: efficient clustering of large molecular libraries.'

setup(
    name='bitbirch',
    version=VERSION,
    description=DESCRIPTION,
    url='https://github.com/mqcomplab/bitbirch.git',  # Update with the actual URL if available
    packages=find_packages(),
    install_requires=[
        'numpy<2,>=1.20',
        'matplotlib',
        'pandas',
        'rdkit',
        'scipy',
        'seaborn',
        'scikit-learn',
        'plotly',
        'iSIM @ git+https://github.com/mqcomplab/iSIM.git',
        'useful_rdkit_utils @ git+https://github.com/PatWalters/useful_rdkit_utils.git'
    ]
)