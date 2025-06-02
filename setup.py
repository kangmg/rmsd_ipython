from setuptools import setup, find_namespace_packages

setup(
    name='rmsd_ipython',
    version='0.0.1',
    author='Kang mingi',
    author_email='kangmg@korea.ac.kr',
    description="Calculate root-mean-square deviation (RMSD) between two sets of XYZ coordinates",
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',  
    url='https://github.com/kangmg/rmsd_ipython',
    keywords=['chemistry','computational chemistry','similarity', 'rmsd'],
    packages=find_namespace_packages(),
    install_requires=[
        "numpy",
        "scipy",
        "matplotlib",
        "pandas",
        "seaborn",
        "ase"
    ],
    classifiers=[ 
        'Programming Language :: Python :: 3',
        "License :: OSI Approved :: BSD License",
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Chemistry'
    ],
    python_requires='>=3.8.0',
)
