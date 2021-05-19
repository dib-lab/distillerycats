from setuptools import setup, find_packages

# read the contents of your README file
from os import path

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()


CLASSIFIERS = [
    "Environment :: Console",
    "Environment :: MacOS X",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: BSD License",
    "Natural Language :: English",
    "Operating System :: POSIX :: Linux",
    "Operating System :: MacOS :: MacOS X",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
    "Programming Language :: Python :: 3.9",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

setup(
    name = 'distillerycats',
    description="Find disease associations across metagenomes with k-mers using sourmash, and then recover pangenomic accessory elements using spacegraphcats.",
    url="https://github.com/dib-lab/distillerycats",
    author="Taylor Reiter, N. Tessa Pierce and C. Titus Brown",
    author_email="tereiter@ucdavis.edu,ntpierce@gmail.com,titus@idyll.org",
    license="BSD 3-clause",
    packages = find_packages(),
    classifiers = CLASSIFIERS,
    entry_points = {'console_scripts': ['distillerycats  = distillerycats.__main__:main']},
    include_package_data=True,
    package_data = { "distillerycats": ["Snakefile", "*.yaml", "*.yml", "*.ipynb"] },
    setup_requires = [ "setuptools>=38.6.0",
                       'setuptools_scm',
                       'setuptools_scm_git_archive',
                        'pytest-runner'],
    use_scm_version = {"write_to": "distillerycats/version.py"},
    install_requires = ['snakemake>=6.3.0', 'click>=8', 'pandas'],
    tests_require=["pytest>=5.1.2", "pytest-dependency"],
    long_description=long_description,
    long_description_content_type="text/markdown",
)
