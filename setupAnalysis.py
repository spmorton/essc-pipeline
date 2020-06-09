#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#==============================================================================
#     setupAnalysis.py
#     Copyright (C) 2020  Scott P Morton (spm3c at mtmail.mtsu.edu)
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
# 
#==============================================================================

"""
Created on Mon Feb 18 18:47:13 2019

@author: Scott P. Morton
Center for Computational Science
Middle Tennessee State University
"""

import os, sys

message = "./setup.py 'the/Analysis/Directory'\n"

analysisDirs = ["Analysis/Additional_Plots", # contains any additional plots such as EVM screeplot
                "Analysis/BE", # Contains binding energies plots and data
                "Analysis/BESI", # Contains BESI related plots
                "Analysis/Datasets", # All correlated data saved in various formats
                "Analysis/EVM_images", # EVM images
                "Analysis/FingerPrints", # Electrophoretic Fingerprints (Stieh et al.)
                "Analysis/HXB2", # HXB2 alignment data
                "Analysis/Logos", # Logos plot data
                "Analysis/Residues", # Residue plots
                "Analysis/seqfiles", # sequence files
                "Analysis/VLoops", # Vloop plots and data
                "Analysis/Types", # Place your types file here
                "Analysis/PCA",
                "Analysis/PCA/FingerPrints"]  

def _cp(src,dest):
    try:
        os.system('cp -r %s %s' % (src,dest))
    except Exception as  err:
        print('(_cp)OS error: {0}'.format(err))
def _Mkdir(path):
    try:
        os.makedirs(path)
    except Exception as  err:
        print('(_Mkdir)OS error: {0}'.format(err))
def _SetWD(path):
    try:
        os.chdir(path)
        print('Current working dir is: %s' % (path))
    except Exception as  err:
        print('(_SetWD)OS error: {0}'.format(err))
def ProcessCMD(cmd,args, message):
    if len(message) > 0:
        print(message)
    wDir = os.getcwd()   
    rc = os.fork()
    if rc == -1:
        print('fork failed')
        raise OSError('os.fork() failed')
    elif rc == 0:
        _SetWD(wDir)
        sys.stdout.flush()
        try:
            os.execvp(cmd,args)
        except Exception as err:
            print('---Child process unable to execute command "%s"\n, Error is-%s' % (cmd,err))
    else:
        status = os.waitpid(rc,0)
        if not(status[1] == 0):
            sys.exit(-1)
            #return status[1]
        else:
            print('Command %s with message %s completed\n' % (cmd,message))
            return 0

def main():
    print(sys.argv[0][0:sys.argv[0].rfind('/')])
    copyExample = False
    root = ''
    if(len(sys.argv) > 1):
        for i in range(len(sys.argv)):
            if(sys.argv[i].upper() == '--HELP' or sys.argv[i].upper() == '-H'):
                print(message)
                sys.exit(0)
            else:
                root = sys.argv[i]
    elif(len(sys.argv) == 1) :
        print(message)
        sys.exit(0)
    _Mkdir(root)
    for folder in analysisDirs:
        _Mkdir('%s/%s' % (root,folder))

        
        
    
if __name__ == "__main__":
    main()      
