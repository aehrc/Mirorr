import os
import sys
from glob import glob
from setuptools import setup, find_packages
from setuptools.extension import Extension
import numpy

setup(
    name='mirorr',
    description='Python tools for Mirorr',
    author='CSIRO AEHRC',
    version='0.0.1',
    license='CSIRO-OSSL',
    packages=find_packages(),
    install_requires=['nipype'],
)
