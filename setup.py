#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 18:23:29 2020

@author: rodrigo
"""

import setuptools
import os
from sys import platform


if platform == "darwin":
    os.environ["CC"] = "gcc-9" 
    os.environ["CXX"] = "g++-9"
        
elif platform == "win32":
    raise ValueError('Windows is not supported yet.')



projectionDDb = setuptools.Extension('projectionDDb',
                                      sources = ['pydbt/sources/projectionDDb/projectionDDb.cpp'],
                                      language='c++',
                                      extra_compile_args = ['-fopenmp'],
                                      extra_link_args=['-lgomp'])

backprojectionDDb = setuptools.Extension('backprojectionDDb',
                                      sources = ['pydbt/sources/backprojectionDDb/backprojectionDDb.cpp'],
                                      language='c++',
                                      extra_compile_args = ['-fopenmp'],
                                      extra_link_args=['-lgomp'])


with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pyDBT",
    version="0.0.1",
    author="Rodrigo Vimieiro",
    description="This package is a python extension of the DBT toolbox from LAVI-USP",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/LAVI-USP/pyDBT",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=["numpy", "matplotlib", "pydicom"],
    ext_modules = [projectionDDb, backprojectionDDb]
)