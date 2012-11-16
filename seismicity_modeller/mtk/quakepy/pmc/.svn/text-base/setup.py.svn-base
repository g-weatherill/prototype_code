# -*- coding: utf-8 -*-
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy

ext_modules=[ 
    Extension( "PMCFunctions", ["PMCFunctions.pyx"] )

# To test this, make libmymath.so or libmymath.dylib from mymath.c
#    Extension("integrate_mymath", ["integrate_mymath.pyx"],
#              include_dirs=['.'],
#              include_libraries=['.'],
#              libraries=['mymath']),

]

setup(
  cmdclass = { 'build_ext': build_ext },
  include_dirs = [ numpy.get_include() ],
  ext_modules = ext_modules,
)
