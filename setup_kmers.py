import numpy
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension(
    name='kmers',
    sources=["count_kmers.pyx"],
    include_dirs=[numpy.get_include()],
    language="c++"
    )]

setup(
    name = 'kmers',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules,
    )
