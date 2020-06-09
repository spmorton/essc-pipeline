# -*- coding: utf-8 -*-
"""
#---------------------------
#
#   Name:libpyzfp.pyx 
#    
#   Cython Script
#   Author: Scott P Morton (spm3c at mtmail.mtsu.edu)
# 
#   libpyzfp is a cython wrapper for ZFP
#
#---------------------------

See ZFP_LICENSE in the libpyzfp source folder for details on ZFP
"""

import cython
import numpy as np
cimport numpy as np

from libc.stdlib cimport malloc, free

cdef extern from "libzfp.h":
    cdef int compress(int argc, char* argv[], void* carray)

def Version():
    return '0.5.6 libpyzfp'

def ZFPC(args,array):
    """ 
        ========
        Compress
        ========
        
        zfp.ZFPC(args,array)
        
        Compress the passed array, per the passed args
        
        Parameters
        ----------        
        args as a string i.e. '-b 1 -z d0.zfp -o d0.rec -d -3 257 161 193 -a 1e-1 -s'
        
        array as a numpy array of data to compress
        
        General options:
        ----------------        
          -h : read/write array and compression parameters from/to compressed header
          
          -q : quiet mode; suppress output
          
          -s : print error statistics
          
        Input and output:
        -----------------        
          -i <path> : uncompressed binary input file (\-\ for stdin)
                     
          -o <path> : decompressed binary output file (\-\ for stdout)
          
          -z <path> : compressed input (w/o -i) or output file (\-\ for stdin/stdout)
          
          -b 0/1    : pass numpy array to/from cython // added by spm
          
        Array type and dimensions (needed with -i):
        -------------------------------------------        
          -f : single precision (float type)
          
          -d : double precision (double type)
          
          -1 <nx> : 1D array dimensions
          
          -2 <nx> <ny> : 2D array dimensions
          
          -3 <nx> <ny> <nz> : 3D array dimensions
          
        Compression parameters (needed with -i):
        ----------------------------------------        
          -r <rate> : fixed rate (# compressed bits per floating-point value)
          
          -p <precision> : fixed precision (# uncompressed bits per value)
          
          -a <tolerance> : fixed accuracy (absolute error tolerance)
          
          -c <minbits> <maxbits> <maxprec> <minexp> : advanced usage
          
              minbits : min # bits per 4^d values in d dimensions
              
              maxbits : max # bits per 4^d values in d dimensions (0 for unlimited)
              
              maxprec : max # bits of precision per value (0 for full)
              
              minexp : min bit plane # coded (-1074 for all bit planes)
                    
        Examples:
        ---------
          >>> -b 1 -z zfile -d -3 257 161 193 -a 1e-1: pass 3D numpy array and commpress to file
          with a max error rate of 0.1
          
          >>> -i file : read uncompressed file and compress to memory
          
          >>> -z file : read compressed file and decompress to memory
          
          >>> -i ifile -z zfile : read uncompressed ifile, write compressed zfile
          
          >>> -z zfile -o ofile : read compressed zfile, write decompressed ofile
          
          >>> -i ifile -o ofile : read ifile, compress, decompress, write ofile
          
          >>> -i file -s : read uncompressed file, compress to memory, print stats
          
          >>> -i - -o - -s : read stdin, compress, decompress, write stdout, print stats
          
          >>> -f -3 100 100 100 -r 16 : 2x fixed-rate compression of 100x100x100 floats
          
          >>> -d -1 1000000 -r 32 : 2x fixed-rate compression of 1M doubles
          
          >>> -d -2 1000 1000 -p 32 : 32-bit precision compression of 1000x1000 doubles
          
          >>> -d -1 1000000 -a 1e-9 : compression of 1M doubles with < 1e-9 max error
          
          >>> -d -1 1000000 -c 64 64 0 -1074 : 4x fixed-rate compression of 1M doubles
    """
    A = array.copy()
    cdef void* cArray = <void*> np.PyArray_DATA(A)

    argsList = args.split(' ')
    cdef char **string_buf = <char **>malloc(len(argsList) * sizeof(char*))
    for i in range(len(argsList)):
        string_buf[i] = argsList[i]
    argcP = len(argsList)

    rc = compress(argcP,string_buf,cArray)

    free(string_buf)

    if(rc):
        print "Action rc = ",rc
        print "Apparently there was an issue, check your args string and parameters\n"
        return -1

    
def ZFPD(args,array):
    """ 
        ========
        DeCompress
        ========
        
        zfp.ZFPD(args,array)
        
        DeCompress and return to the passed array per the passed args
        
        Parameters
        ----------
    
        Decompress
        
        args as a string i.e. '-b 0 -z d0.zfp -o d00.rec -d -3 257 161 193 -a 1e-1'
        
        array as a numpy array as the target of the decompressed data
        
        returns the data recovered as a numpy array


        General options:
        ----------------
        
          -h : read/write array and compression parameters from/to compressed header
          
          -q : quiet mode; suppress output
          
          -s : print error statistics
          
        Input and output:
        -----------------
        
          -i <path> : uncompressed binary input file (\-\ for stdin)
                     
          -o <path> : decompressed binary output file (\-\ for stdout)
          
          -z <path> : compressed input (w/o -i) or output file (\-\ for stdin/stdout)
          
          -b 0/1    : pass numpy array to/from cython // added by spm
          
        Array type and dimensions (needed with -i):
        -------------------------------------------
        
          -f : single precision (float type)
          
          -d : double precision (double type)
          
          -1 <nx> : 1D array dimensions
          
          -2 <nx> <ny> : 2D array dimensions
          
          -3 <nx> <ny> <nz> : 3D array dimensions
          
        Compression parameters (needed with -i):
        ----------------------------------------
          -r <rate> : fixed rate (# compressed bits per floating-point value)
          
          -p <precision> : fixed precision (# uncompressed bits per value)
          
          -a <tolerance> : fixed accuracy (absolute error tolerance)
          
          -c <minbits> <maxbits> <maxprec> <minexp> : advanced usage
          
              minbits : min # bits per 4^d values in d dimensions
              
              maxbits : max # bits per 4^d values in d dimensions (0 for unlimited)
              
              maxprec : max # bits of precision per value (0 for full)
              
              minexp : min bit plane # coded (-1074 for all bit planes)
                    
        Examples:
        ---------
          >>> -b 0 -z zfile -d -3 257 161 193 -a 1e-1: pass 3D numpy array, return decompressed array
          
          >>> -i file : read uncompressed file and compress to memory
          
          >>> -z file : read compressed file and decompress to memory
          
          >>> -i ifile -z zfile : read uncompressed ifile, write compressed zfile
          
          >>> -z zfile -o ofile : read compressed zfile, write decompressed ofile
          
          >>> -i ifile -o ofile : read ifile, compress, decompress, write ofile
          
          >>> -i file -s : read uncompressed file, compress to memory, print stats
          
          >>> -i - -o - -s : read stdin, compress, decompress, write stdout, print stats
          
          >>> -f -3 100 100 100 -r 16 : 2x fixed-rate compression of 100x100x100 floats
          
          >>> -d -1 1000000 -r 32 : 2x fixed-rate compression of 1M doubles
          
          >>> -d -2 1000 1000 -p 32 : 32-bit precision compression of 1000x1000 doubles
          
          >>> -d -1 1000000 -a 1e-9 : compression of 1M doubles with < 1e-9 max error
          
          >>> -d -1 1000000 -c 64 64 0 -1074 : 4x fixed-rate compression of 1M doubles
    """
    A = array.copy()
    cdef void* cArray = <void*> np.PyArray_DATA(A)

    argsList = args.split(' ')
    cdef char **string_buf = <char **>malloc(len(argsList) * sizeof(char*))
    for i in range(len(argsList)):
        string_buf[i] = argsList[i]
    argcP = len(argsList)

    rc = compress(argcP,string_buf,cArray)
    if(rc):
        print "Action rc = ",rc
        print "Apparently there was an issue, check your args string and parameters\n"
        return -1

    free(string_buf)
    return A
        

