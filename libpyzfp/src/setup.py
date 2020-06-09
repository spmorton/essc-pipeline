#!/usr/bin/env python2

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import numpy
import os

os.environ["CFLAGS"] = "-fPIC"

ext = Extension('libpyzfp',
                sources=['libzfp.c','libpyzfp.pyx'],
                include_dirs=[numpy.get_include(),'../inc','../array','inline','template'],
                extra_objects=['../lib/libzfp.a'])
setup(name='libpyzfp', ext_modules = cythonize(ext))

#libraries=['../lib/libzfp'],
#sources=['bitstream.c','decode1f.c','decode1d.c','encode1f.c','encode1d.c','decode2f.c','decode2d.c','encode2f.c','encode2d.c','decode3f.c','decode3d.c','encode3f.c','encode3d.c','zfp.c','zzfp.c','zzfppy.pyx'],
