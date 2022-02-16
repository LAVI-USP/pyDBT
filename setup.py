#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 18:23:29 2020

@author: rodrigo
"""

import setuptools
import os

from sys import platform
from pydbt.functions.utilities2setuptools import locate_cuda, get_custom_build_ext

# Find CUDA (NVCC)
cuda_env = locate_cuda()

# Create custom build ext
custom_build_ext = get_custom_build_ext(cuda_env)


if platform == "darwin":
    os.environ["CC"] = "gcc-11" 
    os.environ["CXX"] = "g++-11"
        
elif platform == "win32":
    raise ValueError('Windows is not supported yet.')



projectionDD = setuptools.Extension('projectionDD',
                                     sources = ['pydbt/sources/projectionDD/projectionDD.cpp'],
                                     language='c++',
                                     extra_compile_args = {"gcc"  : [],
                                                           "nvcc" : []})

backprojectionDD = setuptools.Extension('backprojectionDD',
                                     sources = ['pydbt/sources/backprojectionDD/backprojectionDD.cpp'],
                                     language='c++',
                                     extra_compile_args = {"gcc"  : [],
                                                           "nvcc" : []})                                     


projectionDDb = setuptools.Extension('projectionDDb',
                                      sources = ['pydbt/sources/projectionDDb/projectionDDb.cpp'],
                                      language='c++',
                                      extra_compile_args = {"gcc"  : ['-fopenmp'],
                                                            "nvcc" : []},
                                      extra_link_args=['-lgomp'])

backprojectionDDb = setuptools.Extension('backprojectionDDb',
                                      sources = ['pydbt/sources/backprojectionDDb/backprojectionDDb.cpp'],
                                      language='c++',
                                      extra_compile_args = {"gcc"  : ['-fopenmp'],
                                                            "nvcc" : []},
                                      extra_link_args=['-lgomp'])

projectionDDb_cuda = setuptools.Extension('projectionDDb_cuda',
                                      sources = ['pydbt/sources/projectionDDb_cuda/projectionDDb_cuda.cpp', 'pydbt/sources/projectionDDb_cuda/kernel.cu'],
                                      library_dirs=[cuda_env['lib64']],
                                      libraries=['cudart'],
                                      language='c++',
                                      runtime_library_dirs=[cuda_env['lib64']],
                                      extra_compile_args = {"gcc"  : [],
                                                            "nvcc" : ['-arch=sm_75', '--ptxas-options=-v', '-c', '--compiler-options', "'-fPIC'"]},
                                      include_dirs = [cuda_env['include'], 'cuda-samples/Common/'])

backprojectionDDb_cuda = setuptools.Extension('backprojectionDDb_cuda',
                                      sources = ['pydbt/sources/backprojectionDDb_cuda/backprojectionDDb_cuda.cpp', 'pydbt/sources/backprojectionDDb_cuda/kernel.cu'],
                                      library_dirs=[cuda_env['lib64']],
                                      libraries=['cudart'],
                                      language='c++',
                                      runtime_library_dirs=[cuda_env['lib64']],
                                      extra_compile_args = {"gcc"  : [],
                                                            "nvcc" : ['-arch=sm_75', '--ptxas-options=-v', '-c', '--compiler-options', "'-fPIC'"]},
                                      include_dirs = [cuda_env['include'], 'cuda-samples/Common/'])


with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pyDBT",
    version="0.0.3",
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
    ext_modules = [projectionDD, backprojectionDD, projectionDDb, backprojectionDDb, projectionDDb_cuda, backprojectionDDb_cuda],
    
    # Inject our custom trigger
    cmdclass={'build_ext': custom_build_ext},
)