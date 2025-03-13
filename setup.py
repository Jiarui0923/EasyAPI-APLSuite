from setuptools import setup, find_packages
import os


LONG_DESCRIPTION = '''
This project build RESTful-API based on EasyAPI for APL and its components.
'''
VERSION = '2.0.11'
NAME = 'EasyAPI-APLSuite'


install_requires = [
    "easyapi",
    "biopandas",
    "biopython",
    "freesasa",
    "requests",
    "gpucorex",
    "numpy",
    "pandas",
    "python_on_whales",
    "scipy",
    "Levenshtein",
]


setup(
    name=NAME,
    description="APLSuites for EasyAPI.",
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    version=VERSION,
    packages=find_packages(),
    install_requires=install_requires,
    url="https://git.tulane.edu/apl/easyapi-aplsuite",
    author="Jiarui Li, Marco K. Carbullido, Jai Bansal, Samuel J. Landry, Ramgopal R. Mettu",
    author_email=("jli78@tulane.edu"),
    maintainer=("Jiarui Li"),
    maintainer_email="jli78@tulane.edu",
    license="GPLv3",
    classifiers=[
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX",
        "Operating System :: OS Independent",
        "License :: OSI Approved :: GPLv3 License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.12",
        "Topic :: Software Development :: Build Tools",
    ],
    zip_safe=True,
    platforms=["any"],
)
