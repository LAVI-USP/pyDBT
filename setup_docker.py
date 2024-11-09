from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

CUDA_LIB = '/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/cuda/12.3/targets/x86_64-linux/lib/'
CUDA_INC = '/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/cuda/12.3/include/'
ARCH = '-arch=sm_75'

class CustomBuildExt(build_ext):
    def build_extensions(self):
        self.compiler.src_extensions.append('.cu')

        for ext in self.extensions:
            print(ext)
            if any(src.endswith('.cu') for src in ext.sources):
                ext.extra_compile_args = [ARCH, '--ptxas-options=-v', '-c', '--compiler-options', '-fPIC']
                ext.extra_link_args = ['-shared', '-lcudart'] 
                self.compiler.set_executable("compiler_so", "nvcc")
                self.compiler.set_executable("linker_so", "nvcc")

            elif 'projectionDDb' in ext.name or 'backprojectionDDb' in ext.name:
                ext.extra_compile_args = ['-fopenmp', '-fPIC']
                ext.extra_link_args = ['-shared', '-lgomp', '-fPIC']
                self.compiler.set_executable("compiler_so", "gcc")
                self.compiler.set_executable("linker_so", "gcc")

            else:
                ext.extra_link_args = ['-shared']
                self.compiler.set_executable("compiler_so", "gcc")
                self.compiler.set_executable("linker_so", "gcc")

        super().build_extensions()

projectionDD = Extension(
    'projectionDD',
    sources=['pydbt/sources/projectionDD/projectionDD.cpp'],
    extra_link_args=['-shared'],
    language='c++'
)

backprojectionDD = Extension(
    'backprojectionDD',
    sources=['pydbt/sources/backprojectionDD/backprojectionDD.cpp'],
    extra_link_args=['-shared'],
    language='c++'
)

projectionDDb = Extension(
    'projectionDDb',
    sources=['pydbt/sources/projectionDDb/projectionDDb.cpp'],
    extra_compile_args=['-fopenmp','-fPIC'],
    extra_link_args=['-shared', '-lgomp'],
    language='c++'
)

backprojectionDDb = Extension(
    'backprojectionDDb',
    sources=['pydbt/sources/backprojectionDDb/backprojectionDDb.cpp'],
    extra_compile_args=['-fopenmp','-fPIC'],
    extra_link_args=['-shared', '-lgomp'],
    language='c++'
)

# CUDA-based extensions
projectionDDb_cuda = Extension(
    'projectionDDb_cuda',
    sources=[
        'pydbt/sources/projectionDDb_cuda/projectionDDb_cuda.cpp',
        'pydbt/sources/projectionDDb_cuda/kernel.cu'  
    ],
    library_dirs=[CUDA_LIB],  
    libraries=['cudart'],
    include_dirs=[CUDA_INC, 'cuda-samples/Common/'],
    extra_link_args=['-lcudart'],
    language='c++'
)

backprojectionDDb_cuda = Extension(
    'backprojectionDDb_cuda',
    sources=[
        'pydbt/sources/backprojectionDDb_cuda/backprojectionDDb_cuda.cpp',
        'pydbt/sources/backprojectionDDb_cuda/kernel.cu' 
    ],
    library_dirs=[CUDA_LIB],
    libraries=['cudart'],
    include_dirs=[CUDA_INC, 'cuda-samples/Common/'],
    extra_link_args=['-lcudart'],
    language='c++'
)

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="pyDBT",
    version="0.0.4",
    author="Rodrigo Vimieiro",
    description="This package is a python extension of the DBT toolbox from LAVI-USP",
    long_description=long_description,
    long_description_content_type="text/markdown",
    ext_modules=[projectionDD, backprojectionDD, projectionDDb, backprojectionDDb],
    cmdclass={'build_ext': CustomBuildExt},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=["numpy", "matplotlib", "pydicom"],
)

setup(
    name="pyDBT",
    version="0.0.4",
    author="Rodrigo Vimieiro",
    description="This package is a python extension of the DBT toolbox from LAVI-USP",
    long_description=long_description,
    long_description_content_type="text/markdown",
    ext_modules=[projectionDDb_cuda, backprojectionDDb_cuda],
    cmdclass={'build_ext': CustomBuildExt},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=["numpy", "matplotlib", "pydicom"],
)
