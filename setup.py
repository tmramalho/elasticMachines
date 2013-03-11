'''
Created on Mar 1, 2013

@author: tiago
'''
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

ext_modules = [Extension("csolver", ["csolver.pyx"], include_dirs=[numpy.get_include()], 
						extra_compile_args=["-funroll-loops", "-ftree-vectorize", "-msse2", "-ftree-vectorizer-verbose=5"])]

setup(
  name = 'C ODE solver',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)