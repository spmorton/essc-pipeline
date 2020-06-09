
//==============================================================================
//     libzfp.h - Header file for Cython
//     Copyright (C) 2020  Scott P Morton (spm3c at mtmail.mtsu.edu)
// 
//     This program is free software: you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.
// 
//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
// 
//     You should have received a copy of the GNU General Public License
//     along with this program.  If not, see <http://www.gnu.org/licenses/>.
// 
//==============================================================================


/*
 * Author - Scott P Morton spm3c at mtmail.mtsu.edu
 *
 *
 * Description -    libzfp.c is a modfied version of code originally written by Dr. Peter Lingstrom
 *                  that allows for passing of a void* array for processing via parent cython code
 *                  ZFP is licensed per the license file ZFP_LICENSE in the source code dir.
 *                  Modifications are notated in the source
 *
 *  This is the header file for the cython implementation
 */



#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "zfp.h"
#include "macros.h"

int compress(int argc, char* args[], void* carray);
