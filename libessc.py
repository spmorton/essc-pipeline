#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    ESSC pipeline Library <libessc.py>
    Copyright (C) 2020  Scott P. Morton (spm3c at mtmail.mtsu.edu)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#---------------------------
#
#   Name:libessc.py  
#   Python Script
#   Author: Scott P Morton (spm3c at mtmail.mtsu.edu)
# 
#   libessc is a set of library functions used to model and extract
#   electrostatic binding characteristics of proteins
#
#---------------------------
"""
citation = 'This code is based on and translated from a set of bash and R scripts\n\
 originally written by Joshua L. Phillips PhD'

__title__  = 'libessc'
__author__ = 'Scott P Morton (spm3c at mtmail.mtsu.edu)'
Description = '\
Electro-Static Surface Charge Pipeline Library\n\n\n \
        Copyright (C) 2018  MIDDLE TENNESSEE STATE UNIVERSITY\n \
        Copyright (C) 2018  Scott P. Morton\n \
        \n\n \
         __  __ \n \
        |  \/  |\n \
        | |\/| |\n \
        | |  | |\n \
        |_|  |_|IDDLE\n \
          _____ \n \
         |_   _|\n \
           | |  \n \
           | |  \n \
           |_| ENNESSEE  \n \
          ____ \n \
         / ___| \n \
         \___ \ \n \
          ___) | \n \
         |____/ TATE\n \
          _   _ \n \
         | | | |\n \
         | | | |\n \
         | |_| |\n \
          \___/ NIVERSITY\n \
\
        \n\n \
        MIDDLE TENNESSEE STATE UNIVERSITY\n\n \
        This program comes with ABSOLUTELY NO WARRANTY;\n \
        This is free software, and you are welcome to redistribute it\n \
        under certain conditions; see <http://www.gnu.org/licenses/> for details.\n\n'

import sys, json, os, textwrap, copy, signal, time, gzip, warnings#, cython,imp
import libpyzfp as zfp
import numpy as np
from multiprocessing import Pool
import subprocess as sp
#import subprocess as subp
from glob import glob
import modeller
import modeller.automodel
from modeller import *
import psize

#TODO - Test protocol change to add solvation to complex structure first and then bind the two 
#       structures together, test on a small scale to see if it is worth the extra compute
#TODO - utilize gzip for residue .sas files
#TODO - bring more operation in from the driver, basically the driver should just pass the work unit into the lib
#TODO - Make argument lists to functions consitent in format and method
#TODO - ??? write drivers in C (lose a lot of flexibility in changes for incremental runs)
#TODO - ??? use C data types

def Version():
    return '1.1.6 libessc.py'

def printInfo():
    print ('\n\nThis is the %s\nlibessc - pipeline version: %s\nWritten by: %s \
    \n\n%s\n\n') % (Description,Version(),__author__,citation)

def __init__():
    # global variables
    printInfo()
    print '\nESSC Library loaded and ready to run!'
    global Global
    global Studies
    global abseqs
    global ABmaps
    global myrank
    global complexIN
    global genericIN
    global pinoffset
    global homomodel
    global cgridsubsetList
    Global = {} #Global Settings for all Studies, set by driver
    Studies = {} #All Studies keys in config file, set by Driver
    abseqs = [] #Global list of seq files in AB folder, set by driver
    ABmaps = {} # Global env to ab mappings, set by driver
    myrank = 0
    pinoffset = 0
    complexIN = '' # will be loaded with the complex.in template from src folder
    genericIN = '' # will be loaded with the generic.in template from src folder
    homomodel = '' # will be loaded with the homomodel.py template from src folder
    cgridsubsetList = []
#GROMACS pinoffset requires per node mapping for MPI in HPC envronments that
# have the capacity to run more than on rank per compute node
# value used for modulus must equal number of logical cores per socket

# Completes the intialization of global variables for the lib, after config
# info is loaded
def completeInit():
    global wDir
    global numModels
    global sourceFiles
    global boundTarget
    global unboundTarget
    global boundABTarget
    global unboundABTarget
    global gp120salign
    global ABsalign
    global dockedTarget
    global limitAPBS
    global states
    global pHgridLen
    global studyName
    studyName = Studies['bulk_mode']['studyName']
    wDir = Studies['bulk_mode']['wdir']
    numModels = int(Studies['bulk_mode']['models'])
    sourceFiles = Studies['bulk_mode']['sourceFiles']
    boundTarget = Studies['bulk_mode']['boundTarget']
    unboundTarget = Studies['bulk_mode']['unboundTarget']
    boundABTarget = Studies['bulk_mode']['boundABTarget']
    unboundABTarget = Studies['bulk_mode']['unboundABTarget']
    gp120salign = Studies['bulk_mode']['gp120salign']
    ABsalign = Studies['bulk_mode']['ABsalign']
    dockedTarget = Studies['bulk_mode']['dockedTarget']
    limitAPBS = int(Global['limitAPBS'])
    states = ['bound','unbound']
    pHgridLen = len(cgridsubsetList)

# Residue Maps
Res1 = {'ASX': 'B','CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

Res3 = {'B': 'ASX', 'C': 'CYS', 'D': 'ASP', 'S': 'SER', 'Q': 'GLN', 'K': 'LYS',
     'I': 'ILE', 'P': 'PRO', 'T': 'THR', 'F': 'PHE', 'N': 'ASN', 
     'G': 'GLY', 'H': 'HIS', 'L': 'LEU', 'R': 'ARG', 'W': 'TRP', 
     'A': 'ALA', 'V':'VAL', 'E': 'GLU', 'Y': 'TYR', 'M': 'MET'}

# OS/ File IO Commands always prefixed with '_'
def _cp(src,dest):
    try:
        os.system('cp %s %s' % (src,dest))
    except Exception as  err:
        print('(_cp)OS error: {0}'.format(err))

def _Exists(Name):
    try:
        value = os.path.exists(Name)
    except Exception as  err:
        print('(_Exists)OS error: {0}'.format(err))
    return value

def _Link(src, dest):
    try:
        os.symlink(src,dest)
    except Exception as  err:
        print('(_Link)OS error: {0}'.format(err))

def _Mkdir(path):
    try:
        os.makedirs(path)
    except Exception as  err:
        print('(_Mkdir)OS error: {0}'.format(err))

def _Remove(fileName):
    try:
        os.remove(fileName)
    except Exception as err:
        pass
        #print('(_Remove)OS error: {0}'.format(err))
    
def _Rename(src, dest):
    try:
        os.rename(src,dest)
    except Exception as  err:
        print('(_Rename)OS error: {0}'.format(err))

def _SetEnvVar(var, val):
    try:
        os.putenv(var,val)
    except Exception as  err:
        print('(_setEnvVar)OS error: {0}'.format(err))

def _SetWD(path):
    try:
        os.chdir(path)
        print 'Current working dir is: %s' % (path)
    except Exception as  err:
        print('(_SetWD)OS error: {0}'.format(err))

# Read various files
def _ReadConfig(seq_env_file):
    try:
        with open ( seq_env_file, 'r' ) as seqf:
            seq_env  =  json.load(seqf)
    except (Exception, SystemError, RuntimeError, TypeError, NameError) as  err:
        print('Sequence config '', seq_env_file, '' not found or json syntax error detected.\n')
        print('(_ReadConfig)OS error: {0}'.format(err))
        sys.exit(-1)
    return seq_env

def _ReadFile(filename):
    try:
        with open(filename,'r') as f:
            contents = f.read()
    except Exception as  err:
        print('(_ReadFile)OS error: {0}'.format(err))
        sys.exit(-1)
    return contents

def _ReadSeqFiles(directory):
    try:
        files = glob('%s/*.seq' % (directory))
        files.sort()
    except Exception as  err:
        print('(_ReadSeqFiles)OS error: {0}'.format(err))
        sys.exit(-1)
    return files

def CalcMeansRatios(argsList):
    """Single execution block"""
    apbsDir = argsList[0]
    state = argsList[1]
    numModels = argsList[2]
    #cgridsubsetList = argsList[3]
    seqDir = argsList[3]
    #pHgridLen = len(cgridsubsetList)
    thisDir = '%s/%s' % (apbsDir,state)
    _SetWD(thisDir)
    totalMean = np.zeros(shape=(numModels,pHgridLen))
    totalRatio = np.zeros(shape=(numModels,pHgridLen))
    for j in range(1,numModels + 1):
        strOfJ = str(j).zfill(2)
        for i in range(pHgridLen):
            prefix,conc,pH = cgridsubsetList[i].split()
            if _Exists('%s/%s/%s/skipAPBS' % (thisDir,strOfJ,prefix)):
                continue
            fileName = '%s/%s/%s/%s' % (thisDir,strOfJ,prefix,prefix)
            try:
                area = np.loadtxt('%s.sasa' % (fileName))
                pots = np.loadtxt('%s.pot' % (fileName))
            except Exception as  err:
                print err
            totalMean[j-1][i] = np.sum(pots)/area
            totalRatio[j-1][i] = np.sum(pots)/np.sum(np.abs(pots))
        try:
            np.savetxt('%s/%s/total.mean' % (thisDir,strOfJ),totalMean[j-1])
            np.savetxt('%s/%s/total.ratio' % (thisDir,strOfJ),totalRatio[j-1])
        except Exception as  err:
            print err
    try:
        np.savetxt('%s/total.mean' % (thisDir),totalMean)
        np.savetxt('%s/total.ratio' % (thisDir),totalRatio)
    except Exception as  err:
        print err
    _Link('../apbs/%s/total.mean' % (state),'%s/results/%s-total.mean' % (seqDir,state))            
    _Link('../apbs/%s/total.ratio' % (state),'%s/results/%s-total.ratio' % (seqDir,state))            

def Calc_Means_Ratios(seq,numModels):
    thisDir = '%s/Structures' % (wDir)
    seqDir = '%s/%s' % (thisDir,seq)
    apbsDir = '%s/apbs' % (seqDir)
    """ Threaded execution block"""
    pool = Pool(2)
    #cgridsubsetList = ReadpHFile(cgridsubsetFile)
    argsList = [[apbsDir,'bound',numModels,\
                seqDir],[apbsDir,'unbound',\
                numModels,seqDir]]
    print('Calculating Means and Ratios')
    pool.map(CalcMeansRatios,argsList)
    pool.close()
    pool.join()
    if (seq.find('Complex-')):
        try:
            boundTM = np.loadtxt('%s/results/bound-total.mean' % (seqDir))
            unboundTM = np.loadtxt('%s/results/unbound-total.mean' % (seqDir))
            boundTR = np.loadtxt('%s/results/bound-total.ratio' % (seqDir))
            unboundTR = np.loadtxt('%s/results/unbound-total.ratio' % (seqDir))
            bsubuTM = np.subtract(boundTM,unboundTM)
            np.savetxt('%s/results/bound-unbound-total.mean' % (seqDir),bsubuTM)
            bsubuTR = np.subtract(boundTR,unboundTR)
            np.savetxt('%s/results/bound-unbound-total.ratio' % (seqDir),bsubuTR)
        except Exception as  err:
            print('(Non-Complex)Calc_Means_Ratios error: {0}'.format(err))
    if not (seq.find('Complex-')):
        try:
            structures = seq.split('Complex-')[1]
            Env,AB = structures.split('___')
            cMean = np.loadtxt('%s/results/bound-total.mean' % (seqDir))
            ubEMean = np.loadtxt('%s/%s/results/unbound-total.mean' % (thisDir,Env))
            bEMean = np.loadtxt('%s/%s/results/bound-total.mean' % (thisDir,Env))
            ubABMean = np.loadtxt('%s/%s/results/unbound-total.mean' % (thisDir,AB))
            bABMean = np.loadtxt('%s/%s/results/bound-total.mean' % (thisDir,AB))
            ComplexsubUBEsubUBAB = np.subtract(cMean,ubEMean)
            ComplexsubUBEsubUBAB = np.subtract(ComplexsubUBEsubUBAB,ubABMean)
            ComplexsubBEsubBAB = np.subtract(cMean,bEMean)
            ComplexsubBEsubBAB = np.subtract(ComplexsubBEsubBAB,bABMean)
            np.savetxt('%s/results/Complex-ubA-ubB.mean' % (seqDir),ComplexsubUBEsubUBAB)
            np.savetxt('%s/results/Complex-bA-bB.mean' % (seqDir),ComplexsubBEsubBAB)
            cMeanR = np.loadtxt('%s/results/bound-total.ratio' % (seqDir))
            ubEMeanR = np.loadtxt('%s/%s/results/unbound-total.ratio' % (thisDir,Env))
            bEMeanR = np.loadtxt('%s/%s/results/bound-total.ratio' % (thisDir,Env))
            ubABMeanR = np.loadtxt('%s/%s/results/unbound-total.ratio' % (thisDir,AB))
            bABMeanR = np.loadtxt('%s/%s/results/bound-total.ratio' % (thisDir,AB))
            ComplexsubUBEsubUBABR = np.subtract(cMeanR,ubEMeanR)
            ComplexsubUBEsubUBABR = np.subtract(ComplexsubUBEsubUBABR,ubABMeanR)
            ComplexsubBEsubBABR = np.subtract(cMeanR,bEMeanR)
            ComplexsubBEsubBABR = np.subtract(ComplexsubBEsubBABR,bABMeanR)
            np.savetxt('%s/results/Complex-ubA-ubB.ratio' % (seqDir),ComplexsubUBEsubUBABR)
            np.savetxt('%s/results/Complex-bA-bB.ratio' % (seqDir),ComplexsubBEsubBABR)
        except Exception as  err:
            print('(Complex)Calc_Means_Ratios error: {0}'.format(err))

def CalcResMnR(seq,numModels,state):
    """Single execution block - per seq/state
    Calculates the means of all residues for a sequence state"""
    currentDir = '%s/Structures' % (wDir)
    seqDir = '%s/%s' % (currentDir,seq)
    thisDir = '%s/apbs/%s' % (seqDir,state)
    prefix,conc,pH = cgridsubsetList[0].split()
    baseNames = []
    _SetWD(thisDir)
    if not(_Exists('%s/residues' % (thisDir))):    
        _Mkdir('%s/residues' % (thisDir))
    if not(_Exists('%s/results/residues' % (seqDir))):
        _Mkdir('%s/results/residues' % (seqDir))
    baseFiles = glob('%s/01/%s/residues/*.pot' % (thisDir,prefix))
    for name in baseFiles:
        baseNames.append(name[name.find('-res-')+5:name.find('.pot')].split('-'))
        resID,resName = name[name.find('-res-')+5:name.find('.pot')].split('-')
        totalMean = np.zeros(shape=(numModels,len(cgridsubsetList)))
        for j in range(1,numModels+1):
            strOfJ = str(j).zfill(2)
            for i in range(len(cgridsubsetList)):
                prefix,conc,pH = cgridsubsetList[i].split()
                if _Exists('%s/%s/%s/skipAPBS' % (thisDir,strOfJ,prefix)):
                    continue
                if not(_Exists('%s/%s/residues' % (thisDir,strOfJ))):
                    _Mkdir('%s/%s/residues' % (thisDir,strOfJ))
                thisFile = '%s-%s-res-%s-%s' % (seq,prefix,resID,resName)
                fileName = '%s/%s/%s/residues/%s' % (thisDir,strOfJ,prefix,thisFile)
                try:
                    area = np.loadtxt('%s.sasa' % (fileName))
                    pots = np.loadtxt('%s.pot' % (fileName))
                except Exception as  err:
                    print 'Load area & pots - error is: %s' % (err)
                if(area == 0.0):
                    totalMean[j-1][i] = 0.0
                else:
                    totalMean[j-1][i] = np.sum(pots)/area
            try:
                np.savetxt('%s/%s/residues/%s_%s_total.mean' % (thisDir,strOfJ,resID,resName),totalMean[j-1])
            except Exception as  err:
                print 'save array for this prefix - error is: %s' % (err)
        try:
            np.savetxt('%s/residues/%s_%s_%s_total.mean' % (thisDir,resID,resName,state),totalMean)
        except Exception as  err:
            print 'save array for this model - error is: %s' % (err)
        if not(_Exists('%s/results/residues/%s_%s_%s-total.mean' % (seqDir,resID,resName,state))):
            _Link('../../apbs/%s/residues/%s_%s_%s_total.mean' % (state,resID,resName,state),
                  '%s/results/residues/%s_%s_%s-total.mean' % (seqDir,resID,resName,state))            

def Calc_Binding_Energies(Complex,numModels):
    thisDir = '%s/Structures/%s' % (wDir,Complex)
    _SetWD(thisDir)
    BEb = np.zeros(shape=(numModels,len(cgridsubsetList)))
    BEu = np.zeros(shape=(numModels,len(cgridsubsetList)))
    for i in range(numModels):
        for x in range(len(cgridsubsetList)):
            prefix,conc,pH = cgridsubsetList[x].split()
            structures = Complex.split('Complex-')[1]
            Env,AB = structures.split('___')
            try:
                iFile = gzip.open('%s/apbs/bound/%s/%s/%s-APBS.log.gz' % (thisDir,str(i+1).zfill(2),prefix,prefix),'rb')
                cplexData = iFile.read()
                iFile.close()
            except Exception as err:
                print 'Calc_Binding_energies - error is: %s' % (err)
                raise Exception ('Calc_Binding_energies - error is: %s' % (err))
            Gc = ExtractEnergies(cplexData,'8')
            Geb = ExtractEnergies(cplexData,'4')
            Geub = ExtractEnergies(cplexData,'6')             
            Ga = ExtractEnergies(cplexData,'2')
            BEb[i][x] = Gc - Geb - Ga
            BEu[i][x] = Gc - Geub - Ga
    np.savetxt('%s/results/bound_BEnergies.dat' % (thisDir),BEb)
    np.savetxt('%s/results/unbound_BEnergies.dat' % (thisDir),BEu)

def ExtractEnergies(data,calcID):
    first = data.split('CALCULATION #%s' % (calcID))[1]
    second = first.split('Total electrostatic energy =')[1]
    value = second.split('kJ/mol')[0]
    return float(value)
    
# APBS calculate potentials
def CalcPots(index,seq,state,gridLine):
    """ Single execution block"""
    prefix,conc,pH = gridLine.split()
    strofIndex = str(index+1).zfill(2)
    thisDir = '%s/Structures/%s/apbs/%s/%s/%s' % (wDir,seq,state,strofIndex,prefix)
    tmpDir = Global['tempDir']
    _SetWD(thisDir)
    if _Exists('%s/skipAPBS' % (thisDir)):
        return
    files = glob('*.sas.bin')
    if(len(files) == 1):
        return
    targetPrefix = '%s%s' % (thisDir,thisDir[thisDir.rfind('/'):])
    fileName = '%s.in' % (targetPrefix)
    cmdLine = 'export OMP_NUM_THREADS="%s"; apbs %s' % (Global['OMPThreads'],fileName)
    rvals = subProcessCMD(cmdLine, 'APBS - calculating potentials for %s' % (fileName))
    writeLogs(rvals,'%s-APBS' % (targetPrefix))
    if (rvals[2]):
        print 'APBS failed to process',fileName,'attempting rerun'
        print 'rc = ',rvals[2]
        print 'APBS Error output is: %s' % (rvals[1])
        rvals = subProcessCMD(cmdLine, 'APBS - calculating potentials for %s' % (fileName))
        if (rvals[2]):
            print 'APBS failed to process',fileName
            print 'rc = ',rvals[2]
            print 'APBS Error output is: %s' % (rvals[1])
            raise Exception('APBS failed to process %s' % (fileName))
    if(Global['antibody'] == '1' and seq.find('Complex-') == 0):
        structures = seq.split('Complex-')[1]
        Env,AB = structures.split('___')
        m1_DXFile = '%s/%s_essc_unbound_essc_%s_essc_%s_essc_%s.dx' % (tmpDir,AB,seq,strofIndex,prefix)
        m2_DXFile = '%s/%s_essc_bound_essc_%s_essc_%s_essc_%s.dx' % (tmpDir,Env,seq,strofIndex,prefix)
        m3_DXFile = '%s/%s_essc_unbound_essc_%s_essc_%s_essc_%s.dx' % (tmpDir,Env,seq,strofIndex,prefix)
        m4_DXFile = '%s/%s_essc_bound_essc_%s_essc_%s.dx' % (tmpDir,seq,strofIndex,prefix)
        _Remove(m1_DXFile)
        _Remove(m2_DXFile)
        _Remove(m3_DXFile)
        CalcSASPots(index,seq,state,gridLine,m4_DXFile)
    else:
        DX = '%s/%s_%s_%s_%s.dx' % (tmpDir,seq,state,strofIndex,prefix)
        CalcSASPots(index,seq,state,gridLine,DX)

def calculate(args): # calculates the passed args from MP pool in CalcSASPots
    return args

def CalcRESPots(index,seq,state,gridLine):
    """ Single execution block"""
    prefix,conc,pH = gridLine.split()
    strofIndex = str(index+1).zfill(2)
    thisDir = '%s/Structures/%s/apbs/%s/%s/%s' % (wDir,seq,state,strofIndex,prefix)
    if _Exists('%s/skipAPBS' % (thisDir)):
        return
    rows = 0
    cols = 0
    _SetWD(thisDir)
    protein = thisDir[thisDir.rfind('/')+1:]
    PQRSIZE = '%s.size' % (protein)
    PQRFile = '%s.pqr' % (protein)
    DX = '%s.dx' % (protein)
    Dx = _ReadFile(DX)
    Size = _ReadFile(PQRSIZE)
    z,y,x = Dx[Dx.find('ts ',Dx.find('gridpositions counts'))+3:Dx.find('\n',Dx.find('gridpositions counts'))].split()
    if(_Exists('%s.dx.ZDX' % (protein))):
        DXin = '%s.dx.ZDX' % (protein)
        print 'Calling the Cython monster... be a good little monster and give me my data back'
        mydata = np.zeros((int(x),int(y),int(z)),dtype=np.float64) 
        args = '-b 0 -z %s -d -3 %s %s %s -a %s' % (DXin,x,y,z,Global['ZFPMaxError'])
        mydata = zfp.ZFPD(args,mydata)
        mydata = np.reshape(mydata,(int(x)*int(y)*int(z)))
    else:
        try:
            mydata = Dx[Dx.find('\n',Dx.rfind('data follows\n'))+1:Dx.find('attribute')].split()
        except Exception as  err:
            print err
            sys.exit(-1)
    _Link('../../../../../../tools/VMD_gen-residues.tcl','VMD_gen-residues.tcl')
    cmdLine = 'vmd -e VMD_gen-residues.tcl -args %s %s %s' % (PQRFile, seq, Global['OMPThreads'])
    rvals = subProcessCMD(cmdLine, 'Executing VMD_gen-residues.tcl for %s/%s' % (thisDir,seq))
    if (rvals[2]):
        print ('CalcRESPots - VMD_gen-residues returned error code: %s' % rvals[2])
        writeLogs(rvals, 'VMD_gen-residues')
        raise Exception ('CalcRESPots - VMD_gen-residues returned error code: %s' % rvals[2])
    inVars = ['FINEGRID','FINEDIM','GRIDCENTER']
    ofInterest = ['Num. fine grid pts. = ','Fine grid dims = ','Center = ']
    for i in range(len(inVars)):
        indexL = Size.find(ofInterest[i])
        indexU = Size.find('\n',indexL)
        inVars[i] = Size[indexL + len(ofInterest[i]):indexU].replace('x ','').replace('A','').rstrip().split()
    myFiles = glob('residues/*.sas')
    for thisFile in myFiles:
        SAS = copy.deepcopy(thisFile)
        temp = np.linspace(0,float(inVars[1][0]),float(inVars[0][0]))
        temp = temp - np.mean(temp)
        tempX = temp + float(inVars[2][0])
        temp = np.linspace(0,float(inVars[1][1]),float(inVars[0][1]))
        temp = temp - np.mean(temp)
        tempY = temp + float(inVars[2][1])
        temp = np.linspace(0,float(inVars[1][2]),float(inVars[0][2]))
        temp = temp - np.mean(temp)
        tempZ = temp + float(inVars[2][2])
        mygrid = [list(tempX),list(tempY),list(tempZ)]
        mydiff = [tempX[1]-tempX[0],tempY[1]-tempY[0],tempZ[1]-tempZ[0]]
        try:
            warnings.filterwarnings('error')
            Sas = np.loadtxt(SAS)
            rows,cols = np.shape(Sas)
            mysasInd = np.zeros((rows,cols))
        except Exception as  err:
            name = thisFile[:thisFile.rfind('.')]
            np.savetxt('%s.pot' % (name),np.zeros((1,1)),fmt='%12.10f')
            warnings.filterwarnings('default')
            continue
        warnings.filterwarnings('default')
        for myk in range(3):
            mysasSorted = np.sort(Sas[:,myk].view())
            mysasSortedix = Sas[:,myk].argsort(0)
            mysasPos = np.zeros(rows)
            myi = 0
            for myj in range(rows):
                while (myi < len(mygrid) and mysasSorted[myj] >= mygrid[myk][myi]):
                    myi = myi + 1
                mysasPos[myj] = (myi + 1) - ((mygrid[myk][myi] - mysasSorted[myj])/mydiff[myk])
            mysasInd[mysasSortedix,myk] = mysasPos
        def mytranslate(x,y):
            index = (((x[:,0] - 1) * len(y[1]) * len(y[2])) + ((x[:,1] - 1) * len(y[2])) + x[:,2]) - 1
            return index.astype(int)
        mysasLow = mysasInd - np.floor(mysasInd)
        mysasHigh = np.ceil(mysasInd) - mysasInd
        results = np.zeros((rows,cols))
        try:
            TASKS = [np.asfarray(np.take(mydata,mytranslate(np.column_stack((np.floor(mysasInd[:,0]),\
                                                                                    np.floor(mysasInd[:,1]),\
                                                                                    np.floor(mysasInd[:,2])))\
                                                                                    ,mygrid)))[:,None] * \
                                                    np.column_stack((mysasLow[:,0],mysasLow[:,1],mysasLow[:,2])),\
                    np.asfarray(np.take(mydata,mytranslate(np.column_stack((np.ceil(mysasInd[:,0]),\
                                                                                    np.floor(mysasInd[:,1]),\
                                                                                    np.floor(mysasInd[:,2])))\
                                                                                    ,mygrid)))[:,None] * \
                                                    np.column_stack((mysasHigh[:,0],mysasLow[:,1],mysasLow[:,2])),\
                    np.asfarray(np.take(mydata,mytranslate(np.column_stack((np.ceil(mysasInd[:,0]),\
                                                                                    np.ceil(mysasInd[:,1]),\
                                                                                    np.floor(mysasInd[:,2])))\
                                                                                    ,mygrid)))[:,None] * \
                                                    np.column_stack((mysasHigh[:,0],mysasHigh[:,1],mysasLow[:,2])),\
                    np.asfarray(np.take(mydata,mytranslate(np.column_stack((np.ceil(mysasInd[:,0]),\
                                                                                    np.ceil(mysasInd[:,1]),\
                                                                                    np.ceil(mysasInd[:,2])))\
                                                                                    ,mygrid)))[:,None] * \
                                                    np.column_stack((mysasHigh[:,0],mysasHigh[:,1],mysasHigh[:,2])),\
                    np.asfarray(np.take(mydata,mytranslate(np.column_stack((np.floor(mysasInd[:,0]),\
                                                                                    np.ceil(mysasInd[:,1]),\
                                                                                    np.ceil(mysasInd[:,2])))\
                                                                                    ,mygrid)))[:,None] * \
                                                    np.column_stack((mysasLow[:,0],mysasHigh[:,1],mysasHigh[:,2])),\
                    np.asfarray(np.take(mydata,mytranslate(np.column_stack((np.floor(mysasInd[:,0]),\
                                                                                    np.floor(mysasInd[:,1]),\
                                                                                    np.ceil(mysasInd[:,2])))\
                                                                                    ,mygrid)))[:,None] * \
                                                    np.column_stack((mysasLow[:,0],mysasLow[:,1],mysasHigh[:,2])),\
                    np.asfarray(np.take(mydata,mytranslate(np.column_stack((np.floor(mysasInd[:,0]),\
                                                                                    np.ceil(mysasInd[:,1]),\
                                                                                    np.floor(mysasInd[:,2])))\
                                                                                    ,mygrid)))[:,None] * \
                                                    np.column_stack((mysasLow[:,0],mysasHigh[:,1],mysasLow[:,2])),\
                    np.asfarray(np.take(mydata,mytranslate(np.column_stack((np.ceil(mysasInd[:,0]),\
                                                                                    np.floor(mysasInd[:,1]),\
                                                                                    np.ceil(mysasInd[:,2])))\
                                                                                    ,mygrid)))[:,None] * \
                                                    np.column_stack((mysasHigh[:,0],mysasLow[:,1],mysasHigh[:,2]))]
            pool = Pool(8)
            resultList = pool.imap(calculate,TASKS)
            pool.close()
            pool.join()
            for r in resultList:
                results = np.add(results, r)
        except Exception as  err:
            print argsList
            print 'CalcRESPots (average charge at a point) : {0}'.format(err)
            sys.exit(-1)
        results = np.mean(results / 4,axis=1)
        try:
            name = thisFile[:thisFile.rfind('.')]
            np.savetxt('%s.pot' % (name),results,fmt='%12.10f')
        except Exception as  err:
            print 'CalcRESPots savetxt error: {0}'.format(err)

def CalcSASPots(index,seq,state,gridLine,DX):
    """ Single execution block """
    prefix,conc,pH = gridLine.split()
    strofIndex = str(index+1).zfill(2)
    thisDir = '%s/Structures/%s/apbs/%s/%s/%s' % (wDir,seq,state,strofIndex,prefix)
    rows = 0
    cols = 0
    pHSAS = '%s/prot-%s.pqr.sas' % (thisDir,pH)
    PQRSIZE = '%s%s.size' % (thisDir,thisDir[thisDir.rfind('/'):])
    mydata = procTMPDX([thisDir,DX])
    if(seq.find('Complex-') >= 0):
        dat = DX.split('%s/' % (Global['tempDir']))[1]
        values = dat.split('_essc_')
        if(values[0] == seq):
            SAS = '%s.sas' % (PQRSIZE.split('.size')[0])
            potOUT = '%s.pot' % (PQRSIZE.split('.size')[0])
        else:
           SAS = '%s/Structures/%s/apbs/%s/%s/%s/%s.sas' \
                   % (wDir,values[0],values[1],
                      values[3],values[4].split('.dx')[0],values[4].split('.dx')[0])
           potOUT = '%s_%s_%s.pot' % (values[0],values[1],values[4].split('.dx')[0])
    else:
        SAS = '%s.sas' % (PQRSIZE.split('.size')[0])
        potOUT = '%s.pot' % (PQRSIZE.split('.size')[0])
    Size = _ReadFile(PQRSIZE)
    inVars = ['FINEGRID','FINEDIM','GRIDCENTER']
    ofInterest = ['Num. fine grid pts. = ','Fine grid dims = ','Center = ']
    for i in range(len(inVars)):
        indexL = Size.find(ofInterest[i])
        indexU = Size.find('\n',indexL)
        inVars[i] = Size[indexL + len(ofInterest[i]):indexU].replace('x ','').replace('A','').rstrip().split()
    temp = np.linspace(0,float(inVars[1][0]),float(inVars[0][0]))
    temp = temp - np.mean(temp)
    tempX = temp + float(inVars[2][0])
    temp = np.linspace(0,float(inVars[1][1]),float(inVars[0][1]))
    temp = temp - np.mean(temp)
    tempY = temp + float(inVars[2][1])
    temp = np.linspace(0,float(inVars[1][2]),float(inVars[0][2]))
    temp = temp - np.mean(temp)
    tempZ = temp + float(inVars[2][2])
    mygrid = [list(tempX),list(tempY),list(tempZ)]
    mydiff = [tempX[1]-tempX[0],tempY[1]-tempY[0],tempZ[1]-tempZ[0]]
    try:
        if _Exists(pHSAS):
            Sas = np.loadtxt(pHSAS)
        else:
            Sas = np.fromfile('%s.bin' % (SAS),dtype='double')
            Sas = Sas.reshape(len(Sas)/3,3)
        rows,cols = np.shape(Sas)
        mysasInd = np.zeros((rows,cols))
    except Exception as  err:
        print ('CalcSASPots - load SAS: {0}'.format(err))
        sys.exit(-1)
    for myk in range(3):
        mysasSorted = np.sort(Sas[:,myk].view())
        mysasSortedix = Sas[:,myk].argsort(0)
        mysasPos = np.zeros(rows)
        myi = 0
        for myj in range(rows):
            while (myi < len(mygrid) and mysasSorted[myj] >= mygrid[myk][myi]):
                myi = myi + 1
            mysasPos[myj] = (myi + 1) - ((mygrid[myk][myi] - mysasSorted[myj])/mydiff[myk])
        mysasInd[mysasSortedix,myk] = mysasPos
    def mytranslate(x,y):
        index = (((x[:,0] - 1) * len(y[1]) * len(y[2])) + ((x[:,1] - 1) * len(y[2])) + x[:,2]) - 1
        return index.astype(int)
    mysasLow = mysasInd - np.floor(mysasInd)
    mysasHigh = np.ceil(mysasInd) - mysasInd
    results = np.zeros((rows,cols))
    try:
        TASKS = [np.asfarray(np.take(mydata,mytranslate(np.column_stack((np.floor(mysasInd[:,0]),\
                                                                                np.floor(mysasInd[:,1]),\
                                                                                np.floor(mysasInd[:,2])))\
                                                                                ,mygrid)))[:,None] * \
                                                np.column_stack((mysasLow[:,0],mysasLow[:,1],mysasLow[:,2])),\
                np.asfarray(np.take(mydata,mytranslate(np.column_stack((np.ceil(mysasInd[:,0]),\
                                                                                np.floor(mysasInd[:,1]),\
                                                                                np.floor(mysasInd[:,2])))\
                                                                                ,mygrid)))[:,None] * \
                                                np.column_stack((mysasHigh[:,0],mysasLow[:,1],mysasLow[:,2])),\
                np.asfarray(np.take(mydata,mytranslate(np.column_stack((np.ceil(mysasInd[:,0]),\
                                                                                np.ceil(mysasInd[:,1]),\
                                                                                np.floor(mysasInd[:,2])))\
                                                                                ,mygrid)))[:,None] * \
                                                np.column_stack((mysasHigh[:,0],mysasHigh[:,1],mysasLow[:,2])),\
                np.asfarray(np.take(mydata,mytranslate(np.column_stack((np.ceil(mysasInd[:,0]),\
                                                                                np.ceil(mysasInd[:,1]),\
                                                                                np.ceil(mysasInd[:,2])))\
                                                                                ,mygrid)))[:,None] * \
                                                np.column_stack((mysasHigh[:,0],mysasHigh[:,1],mysasHigh[:,2])),\
                np.asfarray(np.take(mydata,mytranslate(np.column_stack((np.floor(mysasInd[:,0]),\
                                                                                np.ceil(mysasInd[:,1]),\
                                                                                np.ceil(mysasInd[:,2])))\
                                                                                ,mygrid)))[:,None] * \
                                                np.column_stack((mysasLow[:,0],mysasHigh[:,1],mysasHigh[:,2])),\
                np.asfarray(np.take(mydata,mytranslate(np.column_stack((np.floor(mysasInd[:,0]),\
                                                                                np.floor(mysasInd[:,1]),\
                                                                                np.ceil(mysasInd[:,2])))\
                                                                                ,mygrid)))[:,None] * \
                                                np.column_stack((mysasLow[:,0],mysasLow[:,1],mysasHigh[:,2])),\
                np.asfarray(np.take(mydata,mytranslate(np.column_stack((np.floor(mysasInd[:,0]),\
                                                                                np.ceil(mysasInd[:,1]),\
                                                                                np.floor(mysasInd[:,2])))\
                                                                                ,mygrid)))[:,None] * \
                                                np.column_stack((mysasLow[:,0],mysasHigh[:,1],mysasLow[:,2])),\
                np.asfarray(np.take(mydata,mytranslate(np.column_stack((np.ceil(mysasInd[:,0]),\
                                                                                np.floor(mysasInd[:,1]),\
                                                                                np.ceil(mysasInd[:,2])))\
                                                                                ,mygrid)))[:,None] * \
                                                np.column_stack((mysasHigh[:,0],mysasLow[:,1],mysasHigh[:,2]))]
        pool = Pool(4)        
        resultList = pool.imap(calculate,TASKS)
        pool.close()
        pool.join()
        for r in resultList:
            results = np.add(results, r)
    except Exception as  err:
        print argsList
        print 'CalcSASPots (average charge at a point) : {0}'.format(err)
        sys.exit(-1)
    results = np.mean(results / 4,axis=1)
    try:    
        np.savetxt(potOUT,results,fmt='%12.10f')
    except Exception as  err:
        print 'CalcSASPots savetxt error: {0}'.format(err)
        sys.exit(-1)
    if not(_Exists('%s.bin' % (SAS))):
        try:
            Sas.tofile('%s.bin' % (SAS))
            _Remove(pHSAS)
        except Exception as  err:
            print 'CalcSASPots binary to file error: {0}'.format(err)
            sys.exit(-1)
    myFiles = ['io.mc','prot-%s.propka' % (pH),]
    for thisFile in myFiles:
        if (_Exists(thisFile)):
            _Remove(thisFile)

# Convert raw sequence into Modeller format
def Convert_Seq(src, dest):
    """ converts the sequence file to modeller format for processing"""
    try:
        seq = open( src, 'r' )
        sequence = seq.read()
        seq.close()
    except Exception as  err:
        print("Sequence config '%s' not found.\n" % (src))
        print('OS error: {0}'.format(err))
    try:
        sequence = sequence.replace('\n','')
        modellerFile = open (dest, 'w')
        modellerFile.write('>P1;model\nsequence:model:1:A:%s:A:::0.00: 0.00\n' % (str(len(sequence))))
        modellerFile.write(textwrap.fill(sequence, 50))
        modellerFile.write('*\n')
        modellerFile.close()
    except Exception as  err:
        print('\nFile creation error for modellerFile\n\n')
        print('OS error: {0}'.format(err))

def CreateDirectoryStructure(seq):
    print 'Creating directory structures for %s' % (seq)
    thisDir = '%s/Structures/%s' % (wDir,seq)
    _Mkdir('%s/plots' % (thisDir))
    _Mkdir('%s/results' % (thisDir))
    _Mkdir('%s/modeller' % (thisDir))
    for index in range(1,  numModels + 1):
        _Mkdir('%s/frodan/%s/bound' % (thisDir,str(index).zfill(2)))
        _Mkdir('%s/frodan/%s/unbound' % (thisDir,str(index).zfill(2)))
        for j in range(len(cgridsubsetList)):
            _Mkdir('%s/apbs/bound/%s/%s' % (thisDir,str(index).zfill(2),cgridsubsetList[j].split()[0]))
            _Mkdir('%s/apbs/unbound/%s/%s' % (thisDir,str(index).zfill(2),cgridsubsetList[j].split()[0]))
def R_Analysis():
    _SetWD(wDir)
    cmdLine = 'Rscript %s/BESI_EVM_analysis.r "%s" %s' % (wDir, studyName, wDir)
    rvals = subProcessCMD(cmdLine, 'BESI/EVM analysis')
    if (rvals[2]):
        writeLogs(rvals,'BESI_EVM')
        print ('BESI/EVM analysis generated and error')
        raise Exception ('BESI/EVM returned error code: %s' % (rvals[2]))
def genComplexIN(index,Complex,gridLine,PDB='prot.pdb'):
    strofIndex = str(index+1).zfill(2)
    prefix,conc,pH = gridLine.split(' ')
    structures = Complex.split('Complex-')[1]
    Env,AB = structures.split('___')
    IN = '%s.in' % (prefix)
    DX = prefix
    tmpDir = Global['tempDir']
    PQR = '%s-%s.pqr' % (PDB.split('.')[0],pH)
    PQRSIZE = '%s.size' % (PQR)
    PCONC = format(float(conc), '.6f')
    NCONC = PCONC
    mol1 = '../../../../../%s/apbs/unbound/%s/%s/%s.pqr' % (AB,strofIndex,prefix,prefix)
    mol2 = '../../../../../%s/apbs/bound/%s/%s/%s.pqr' % (Env,strofIndex,prefix,prefix)
    mol3 = '../../../../../%s/apbs/unbound/%s/%s/%s.pqr' % (Env,strofIndex,prefix,prefix)
    mol4 = '%s-%s.pqr' % (PDB.split('.')[0],pH)
    m1_DXFile = '%s/%s_essc_unbound_essc_%s_essc_%s_essc_%s' % (tmpDir,AB,Complex,strofIndex,DX)
    m2_DXFile = '%s/%s_essc_bound_essc_%s_essc_%s_essc_%s' % (tmpDir,Env,Complex,strofIndex,DX)
    m3_DXFile = '%s/%s_essc_unbound_essc_%s_essc_%s_essc_%s' % (tmpDir,Env,Complex,strofIndex,DX)
    m4_DXFile = '%s/%s_essc_bound_essc_%s_essc_%s' % (tmpDir,Complex,strofIndex,DX)
    _SetWD('%s/Structures/%s/apbs/bound/%s/%s' % (wDir,Complex,strofIndex,prefix))
    if not(_Exists(IN)):
        info = _ReadFile(PQRSIZE)
        inVars = ['GRID','COARSEGRID','FINEGRID','GRIDCENTER']
        ofInterest = ['Num. fine grid pts. = ','Coarse grid dims = ', \
                        'Fine grid dims = ','Center = ']
        for i in range(len(inVars)):
            indexL = info.find(ofInterest[i])
            indexU = info.find('\n',indexL)
            inVars[i] = info[indexL + len(ofInterest[i]):indexU].replace('x ','').replace('A','').rstrip()
        outFile = complexIN.replace('++!!file1!!++',mol1)
        outFile = outFile.replace('++!!file2!!++',mol2)
        outFile = outFile.replace('++!!file3!!++',mol3)
        outFile = outFile.replace('++!!file4!!++',mol4)
        outFile = outFile.replace('++!!var0!!++',inVars[0])
        outFile = outFile.replace('++!!var1!!++',inVars[1])
        outFile = outFile.replace('++!!var2!!++',inVars[2])
        outFile = outFile.replace('++!!var3!!++',inVars[3])
        outFile = outFile.replace('++!!PCONC!!++',PCONC)
        outFile = outFile.replace('++!!NCONC!!++',NCONC)
        outFile = outFile.replace('++!!m1_DXout!!++',m1_DXFile)
        outFile = outFile.replace('++!!m2_DXout!!++',m2_DXFile)
        outFile = outFile.replace('++!!m3_DXout!!++',m3_DXFile)
        outFile = outFile.replace('++!!m4_DXout!!++',m4_DXFile)
        try:
            with open(IN,'w') as f:
                f.write(outFile)
        except Exception as  err:
            print err
            sys.exit(-1)

def genIN(index,seq,state,gridLine, PDB='prot.pdb'):
    """ Reads parameters from psize output and\n
  		writes .in file to target directory\n
  		Generates all needed file links\n"""
    strofIndex = str(index+1).zfill(2)
    apbsDir = '%s/Structures/%s/apbs/%s/%s' % (wDir,seq,state,strofIndex)
    prefix,conc,pH = gridLine.split(' ')
    _SetWD('%s/%s' % (apbsDir,prefix))
    if _Exists('../../../../frodan/%s/%s/skipAPBS' % (strofIndex,state)):
        _Link('../../../../../../src/skipAPBS', 'skipAPBS')        
        return
    IN = '%s.in' % (prefix)
    DX = prefix
    PCONC = format(float(conc), '.6f')
    NCONC = PCONC
    PQR = '%s-%s.pqr' % (PDB.split('.')[0],pH)
    PQRSIZE = '%s.size' % (PQR)
    DXFile = '%s/%s_%s_%s_%s' % (Global['tempDir'],seq,state,strofIndex,DX)
    if not(_Exists(IN)):
        info = _ReadFile(PQRSIZE)
        inVars = ['GRID','COARSEGRID','FINEGRID','GRIDCENTER']
        ofInterest = ['Num. fine grid pts. = ','Coarse grid dims = ', \
                        'Fine grid dims = ','Center = ']
        for i in range(len(inVars)):
            indexL = info.find(ofInterest[i])
            indexU = info.find('\n',indexL)
            inVars[i] = info[indexL + len(ofInterest[i]):indexU].replace('x ','').replace('A','').rstrip()
        outFile = genericIN.replace('++!!file1!!++',PQR)
        outFile = outFile.replace('++!!var0!!++',inVars[0])
        outFile = outFile.replace('++!!var1!!++',inVars[1])
        outFile = outFile.replace('++!!var2!!++',inVars[2])
        outFile = outFile.replace('++!!var3!!++',inVars[3])
        outFile = outFile.replace('++!!PCONC!!++',PCONC)
        outFile = outFile.replace('++!!NCONC!!++',NCONC)
        outFile = outFile.replace('++!!m1_DXout!!++',DXFile)
        try:    
            with open(IN,'w') as f:
                f.write(outFile)
        except Exception as  err:
            print err
            sys.exit(-1)
    
def GenSAS(index,seq,state,gridLine):
    """ Single execution block"""
    prefix,conc,pH = gridLine.split()
    strofIndex = str(index+1).zfill(2)
    thisDir = '%s/Structures/%s/apbs/%s/%s/%s' % (wDir,seq,state,strofIndex,prefix)
    _SetWD(thisDir)
    pqrFile = 'prot-%s.pqr' % (pH)
    if (_Exists('%s/skipAPBS' % (thisDir)) or 
            _Exists('prot-%s.pqr.sas' % (pH)) or 
                _Exists('%s.sas.bin' % (prefix))):
        return
    _Link('../../../../../../tools/VMD_gen-SAS.tcl','VMD_gen-SAS.tcl')
    _Link('prot-%s.pqr.sasa' % (pH), '%s.sasa' % (prefix))
    _Link('prot-%s.pqr.sas' % (pH), '%s.sas' % (prefix))

    cmdLine = 'vmd -e VMD_gen-SAS.tcl -args %s %s' % (pqrFile, Global['OMPThreads'])
    rvals = subProcessCMD(cmdLine, 'Executing VMD_gen-SAS.tcl for %s/%s' % (thisDir,pqrFile))
    if (rvals[2]):
        print ('GenSAS - VMD_gen-SAS returned error code: %s' % (rvals[2]))
        raise Exception ('GenSAS - VMD_gen-SAS returned error code: %s' % (rvals[2]))

def GetLanguage():
    """ Returns the lexica for config json file """
    lexica = {
        "Global":{
            "values":['OMPThreads','asyncThreads','SeqThreads','CompressData',\
                        'ZFPMaxError','mdrunThreads','antibody','ABmode',\
                        'limitAPBS','multiChain','socketsPerNode',\
                        'physicalCores','SMT-HT','GMX-Margin','doRAnalysis'],
            "ZFPError":['1e-1','1e-2','1e-3','1e-4','1e-5','1e-6','1e-7','1e-8','1e-9'],

            },
        "Studies":{
            "dirs": ['wdir','templates','seqFiles','sourceFiles',
                     'ABtemplates','ABseqFiles'],
            "files": ['boundTarget','unboundTarget','cgridsubsetFile',
                      'boundABTarget','unboundABTarget','gp120salign',
                      'ABsalign','dockedTarget'],
            "values":['models']
        }
    }
    return lexica

#TODO FIX THIS FOR NEW METHOD
def Get_New_Model(recList,tag):
    index = recList[0]
    frodanDir = recList[1]
    boundTarget = recList[2]
    unboundTarget = recList[3]
    seq = frodanDir[frodanDir[:frodanDir.rfind('/')].rfind('/')+1:frodanDir.rfind('/')]
    modellerDir = '%sMR_%s_%s' % (frodanDir[:frodanDir.rfind('/')+1],str(index),seq)
    print 'modeller dir is ', modellerDir
    if not(_Exists(modellerDir)):
        _Mkdir(modellerDir)
        print "Model Recovery - Created recovery model folder"
        print "Model Recovery - Recreating model"
        Run_Modeller(seq,modellerDir,'1',tag)
        _cp(modellerDir+"/model.B9999"+str(1).zfill(4)+".pdb","../modeller/model.B9999"+str(index).zfill(4)+".pdb")
        # clean up the mess
        print "Model Recovery - Cleaning up from previous run"
        _SetWD(frodanDir+"/"+str(index).zfill(2))
        for fileName in ['topol.tpr','model.pdb','topol.top',\
                            'conf.gro','grompp.mdp','mdout.mdp',\
                            'posre.itp','confout.gro','conf.box.gro',\
                            'posre_Protein_chain_A.itp','posre_Protein_chain_B.itp',\
                            'topol_Protein_chain_A.itp','topol_Protein_chain_B.itp',\
                            'model.bak.pdb','md.log']:
            _Remove(fileName)
        #reprep frodan
        print "Model Recovery - Running PrepFrodaN after new model built"
        modellerDir = frodanDir[:frodanDir.rfind('/')+1]+'modeller'
        src = frodanDir[:frodanDir.rfind('/')+1]+'src'
        args = []
        args.append(index)
        args.append(frodanDir)
        args.append(modellerDir)
        args.append(src)
        args.append(boundTarget)
        args.append(unboundTarget)
        args.append(boundTarget)
        args.append(unboundTarget)
        PrepFrodaN(args,tag)
        of = open(frodanDir+"/"+str(index).zfill(2)+'/model_replaced','wr')
        #of.truncate()
        of.write('done')
        of.close()
        #_SetWD(currentWDir)
        print "Model Recovery - finished recovery returning"
    # I am the other process, either bound or unbound for this model and we ran asynch
    else:
        print "Model Recovery -waiting for new model from other process"
        while (True):
            time.sleep(30)
            if _Exists(frodanDir+"/"+str(index).zfill(2)+'/model_replaced'):
                print "Model Recovery - finished waiting"
                break

def GenPsize(index,seq,state,gridLine,PDB='prot.pdb'):
    """Generate size information for apbs using psize.py"""
    strofIndex = str(index+1).zfill(2)
    apbsDir = '%s/Structures/%s/apbs/%s/%s' % (wDir,seq,state,strofIndex)
    #'%s/Structures/%s/frodan/%s' % (wDir,seq,str(model+1).zfill(2)),
    prefix,conc,pH = gridLine.split(' ')
    PQR = '%s-%s.pqr' % (PDB.split('.pdb')[0],pH)
    PQRSIZE = '%s.size' % (PQR)
    _SetWD('%s/%s' % (apbsDir,prefix))
    if _Exists('../../../../frodan/%s/%s/skipAPBS' % (strofIndex,state)):
        print('skipped - skipAPBS exists')
        return
    _Link('../../../../frodan/%s/%s/%s.pdb' % (strofIndex,state,state), 'prot.pdb')
    _Link('prot-%s.pqr' % (pH), '%s.pqr' % (prefix))
    _Link('prot-%s.pqr.size' % (pH), '%s.size' % (prefix))
    if not(_Exists(PQR)):
        cmdLine = 'pdb2pqr.py --ff amber --ph-calc-method propka --with-ph=%s %s %s' % (pH,PDB,PQR)
        rvals = subProcessCMD(cmdLine,'Executing pdb2pqr (GenPsize)')
        if (rvals[2]):
            print 'pdb2pqr failed for : %s %s attempting to rerun process' % (PDB,PQR)
            rvals = subProcessCMD(cmdLine,'Executing pdb2pqr (GenPsize)')
            if (rvals[2]):
                writeLogs(rvals, 'pdb2pqr')
                raise Exception('pdb2pqr failed in GenPsize')
    if not(_Exists(PQRSIZE)):
        if(apbsDir.find('Complex-') == -1):
            try:
                outfile = open(PQRSIZE,'wr')
                psz=psize.Psize()
                psz.runPsize(PQR)
                outfile.write('# Constants used: \n');
                for key in psz.constants.keys():
                        outfile.write('# \t%s: %s\n' % (key, psz.constants[key]))
                outfile.write(psz.printResults())
                outfile.flush()
                os.fsync(outfile)
                outfile.close()
            except Exception as err:
                print 'psize error - :', err
        else:
            _Link('../modelPSIZE/merged.psize',PQRSIZE)

def prep_gmx(srcFile,wDir):
    """Prepares PDB for gmx format via pdb2gmx"""
    _SetWD(wDir)
    cmdLine = 'gmx pdb2gmx -ignh -ff amber99sb-ildn -water tip3p -f %s -o conf.gro -ter -merge all' % (srcFile)
    rvals = subProcessCMD(cmdLine,'Executing pdb2gmx for %s' % (srcFile))
    if (rvals[2]):
        print 'pdb2gmx failed for %s attempting rerun' % (srcFile)
        print 'rc is: %s' % (rvals[2])
        print 'stderr has: %s' % (rvals[1])
        rvals = subProcessCMD(cmdLine,'Executing pdb2gmx for %s' % (srcFile))
        if (rvals[2]):
            print 'pdb2gmx failed for %s' % (srcFile)
            writeLogs(rvals,'pdb2gmx')
            raise Exception('Prep frodan - pdb2gmx failed to process')

def gmx_minimize(srcFile,destFile,wDir):
    """performs a gromacs minimization on the source file to the output file"""
    _SetWD(wDir)
    data = ['editconf','editconf']
    cmdLine = 'gmx editconf -f %s -o conf.box.gro -bt triclinic -d %d' % (srcFile,int(Global['GMX-Margin']))
    rvals = subProcessCMD(cmdLine, '(gmx_minimize)Executing editconfig (pre-mdrun)')
    if (rvals[2]):
        print 'editconf failed to process in %s attempting rerun' % (wDir)
        print 'rc is: %s' % (rvals[2])
        rvals = subProcessCMD(cmdLine, '(gmx_minimize)Executing editconfig (pre-mdrun)')
        if (rvals[2]):
            print ' editconf failed to process in %s' % (wDir)
            writeLogs(rvals,'editconf')
            raise Exception('editconf failed to process')
    data[0]  = '%s\n%s' % (data[0],rvals[0])
    data[1]  = '%s\n%s' % (data[1],rvals[1])
    cmdLine = 'gmx grompp -c conf.box.gro -maxwarn 10'
    rvals = subProcessCMD(cmdLine,'(gmx_minimize)Executing grompp')
    if (rvals[2]):
        print 'grompp failed to process in %s attempting rerun' % (wDir)
        print 'rc is: %s' % (rvals[2])
        rvals = subProcessCMD(cmdLine,'(gmx_minimize)Executing grompp')
        if (rvals[2]):
            print 'grompp failed to process in %s' % (wDir)
            writeLogs(rvals,'grompp')
            raise Exception('grompp failed to process')
    data[0]  = '%s\n\n%s\n%s' % (data[0],'grompp',rvals[0])
    data[1]  = '%s\n\n%s\n%s' % (data[1],'grompp',rvals[1])
    pinOffset = ((myrank-1)%int(Global['socketsPerNode']))*int(Global['physicalCores'])
    cmdLine = 'gmx mdrun -nt %s -pin on -pinoffset %s' % (Global['mdrunThreads'],str(pinOffset))
    rvals = subProcessCMD(cmdLine,'(gmx_minimize)Executing mdrun')
    if (rvals[2]):
        print 'mdrun failed to process in %s attempting rerun' % (wDir)
        print 'rc is: %s' % (rvals[2])
        cmdLine = 'gmx mdrun -nt %s -pin on -pinoffset %s' % (Global['mdrunThreads'],str(pinOffset))
        rvals = subProcessCMD(cmdLine,'(gmx_minimize)Executing mdrun')
        if (rvals[2]):
            print 'mdrun failed to process in %s' % (wDir)
            writeLogs(rvals,'mdrun')
            raise Exception('mdrun failed to process')
    data[0]  = '%s\n\n%s\n%s' % (data[0],'mdrun',rvals[0])
    data[1]  = '%s\n\n%s\n%s' % (data[1],'mdrun',rvals[1])
    cmdLine = 'gmx editconf -f confout.gro -label G -o %s' % (destFile)
    rvals = subProcessCMD(cmdLine,'(gmx_minimize)Executing editconf2')
    if (rvals[2]):
        print 'editconf(2) failed to process in %s attmepting rerun' % (wDir)
        print 'rc is: %s' % (rvals[2])
        cmdLine = 'gmx editconf -f confout.gro -label G -o %s' % (destFile)
        rvals = subProcessCMD(cmdLine,'(gmx_minimize)Executing editconf2')
        if (rvals[2]):
            print 'editconf failed to process in %s' % (wDir)
            writeLogs(rvals,'editconf2')
            raise Exception('editconf failed to process')
    data[0]  = '%s\n\n%s\n%s' % (data[0],'editconf2',rvals[0])
    data[1]  = '%s\n\n%s\n%s' % (data[1],'editconf2',rvals[1])
    writeLogs(data,'gmx_minimize')

def handler(sigNum, frame):
    print 'signal handler called with signal', sigNum
    
def Handler_Set():
    for i in range(1,30):
        if i == 9 or i == 19 or i == 29:
            continue
        signal.signal(i,handler)

def multifrodan(model, target):
    num = 0
    prefix = 'frodan%s' % (str(num).zfill(2))
    current = '%s.pdb' %(prefix)
    _cp(model, current)
    num += 1
    # simulation
    data = ['frodanpp','frodanpp']
    cmdLine = 'frodanpp.py -i %s -t %s -o options.xml' % (current,target)
    rvals = subProcessCMD(cmdLine,'Executing frodanpp for %s' % (current))
    if (rvals[2]):
        print 'frodan pre-processor returned error code: %s\n attempting rerun' % (rvals[2])
        rvals = subProcessCMD(cmdLine,'Executing frodanpp for %s' % (current))
        if (rvals[2]):
            writeLogs(rvals,'frodanpp')
            return -1
    data[0] = '%s\n%s' % (data[0],rvals[0])
    data[1] = '%s\n%s' % (data[1],rvals[1])
    data[0] = '%s\n\n%s' % (data[0], 'Please check the following alignment if problems arise...\n PDBaligned.seq\n')
    cmdLine = 'frodan options.xml'
    rvals = subProcessCMD(cmdLine, 'Executing frodan for %s' % (current))
    if (rvals[2]):
        print 'frodan returned error code: %s\n attempting rerun' % (rvals[2])
        rvals = subProcessCMD(cmdLine, 'Executing frodan for %s' % (current))
        if (rvals[2]):
            writeLogs(rvals, 'frodan')
            return -1
    data[0] = '%s\n\n%s\n%s' % (data[0],'frodan',rvals[0])
    data[1] = '%s\n\n%s\n%s' % (data[1],'frodan',rvals[1])
    # trajectory processing
    last = glob('processed_snapshot_*.pdb')
    last.sort()
    cmdLine = 'gmx trjcat -f processed_snapshot_*.pdb -o temp.xtc -cat'
    rvals = subProcessCMD(cmdLine,'Executing trjcat (multifrodan)')
    if (rvals[2]):
        print 'trjcat returned error code: %s\n attempting rerun' % (rvals[2])
        rvals = subProcessCMD(cmdLine,'Executing trjcat (multifrodan)')
        if (rvals[2]):
            writeLogs(rvals,'trjcat')
            return -1
    data[0] = '%s\n\n%s\n%s' % (data[0],'trjcat',rvals[0])
    data[1] = '%s\n\n%s\n%s' % (data[1],'trjcat',rvals[1])
    cmdLine = 'gmx trjconv -f temp.xtc -o %s.xtc -timestep 1' % (prefix) 
    rvals = subProcessCMD(cmdLine,'Executing trjconv')
    if (rvals[2]):
        print 'trjconv returned arror code: %s\n attempting rerun' % (rvals[2])
        rvals = subProcessCMD(cmdLine,'Executing trjconv')
        if (rvals[2]):
            writeLogs(rvals,'trjconv')
            return -1
    data[0] = '%s\n\n%s\n%s' % (data[0],'trjconv',rvals[0])
    data[1] = '%s\n\n%s\n%s' % (data[1],'trjconv',rvals[1])
    _Remove('temp.xtc')
    _cp(last.pop(), 'frodan%s.pdb' % (str(num).zfill(2)))
    #cleanup
    try:
        last = glob('processed_snapshot_*.pdb')
    except Exception as  err:
        print 'Error processing glob', err
    last = last + ['atomtypes.txt', 'cov.txt', 'rc.txt', 'rmsdToTarget.txt', \
                    'targ_atomtypes.txt', 'targcov.txt', 'targmap.txt', \
                    'targrc.txt', 'processed.pdb', 'targprocessed.pdb',\
                    'PDBaligned.seq']
    for item in last:
        if (_Exists(item)):
            _Remove(item)
    cmdLine = 'gmx trjcat -cat -f frodan??.xtc -o temp.xtc'
    rvals = subProcessCMD(cmdLine,'Executing trjcat2')
    if (rvals[2]):
        print 'trjcat returned arror code: %s\n attempting rerun' % (rvals[2])
        rvals = subProcessCMD(cmdLine,'Executing trjcat')
        if (rvals[2]):
            writeLogs(rvals,'trjcat2')
            return -1
    data[0] = '%s\n\n%s\n%s' % (data[0],'trjcat2',rvals[0])
    data[1] = '%s\n\n%s\n%s' % (data[1],'trjcat2',rvals[1])
    cmdLine = 'echo -e "Protein\nSystem" | \
                gmx trjconv -s frodan00.pdb -f temp.xtc -o frodan.xtc\
                -fit rot+trans -timestep 1'
    rvals = subProcessCMD(cmdLine,'Executing trjconv2 (Protein-System)')
    if (rvals[2]):
        print 'trjconv2 returned arror code: %s\n attempting rerun' % (rvals[2])
        rvals = subProcessCMD(cmdLine,'Executing trjconv2')
        if (rvals[2]):
            writeLogs(rvals,'trjconv2')
            return -1
    data[0] = '%s\n\n%s\n%s' % (data[0],'trjconv2',rvals[0])
    data[1] = '%s\n\n%s\n%s' % (data[1],'trjconv2',rvals[1])
    writeLogs(data,'multifrodan')
    _Remove('temp.xtc')
    return 0

def MergeAllPDB(index,Complex, PDB='prot.pdb'):
    structures = Complex.split('Complex-')[1]
    Env,AB = structures.split('___')
    thisDir = '%s/Structures/%s/apbs/bound/%s' % (wDir,Complex,index)
    mol1 = '%s/Structures/%s/frodan/%s/unbound/unbound.pdb' % (wDir,AB,index)
    mol2 = '%s/Structures/%s/frodan/%s/bound/bound.pdb' % (wDir,Env,index)
    mol3 = '%s/Structures/%s/frodan/%s/unbound/unbound.pdb' % (wDir,Env,index)
    mol4 = '%s/Structures/%s/frodan/%s/bound/bound.pdb' % (wDir,Complex,index)
    if not(_Exists('%s/modelPSIZE' % (thisDir))):
        _Mkdir('%s/modelPSIZE' % (thisDir))
    _SetWD('%s/modelPSIZE' % (thisDir))
    if _Exists('merged.psize'):
        return
    _Link('../../../../../../tools/VMD_merge_all_Complex.tcl','mergeAll.tcl')
    cmdLine = 'vmd -e mergeAll.tcl -args %s %s %s %s' % (mol1, mol2, mol3, mol4)
    rvals = subProcessCMD(cmdLine, 'Executing VMD_merge_all_Complex.tcl')
    if  (rvals[2]):
        print 'MergeAllPDB - VMD returned error code: %s' % (rvals[2])
        writeLogs(rvals,'VMDmergeall')
        raise Exception ('MergeAllPDB - VMD returned error code: %s' % (rvals[2]))
    try:
        psz=psize.Psize()
        psz.runPsize('merged.pdb')
        outfile = open('merged.psize','wr')
        outfile.write('# Constants used: \n');
        for key in psz.constants.keys():
                outfile.write('# \t%s: %s\n' % (key, psz.constants[key]))
        outfile.write(psz.printResults())
        outfile.flush()
        os.fsync(outfile)
        outfile.close()
    except Exception as err:
        print 'psize error - :', err

def PrepComplex(index,Complex):
    """Prep the Env and AB into a complex PDB"""
    structures = Complex.split('Complex-')[1]
    Env,AB = structures.split('___')
    strOfIndex = str(index).zfill(2)
    _SetWD('%s/Structures/%s/frodan/%s/unbound/' % (wDir,Complex,strOfIndex))
    if(_Exists('merged.pdb')):
        return
    _Link('../../../../../src/skipAPBS', 'skipAPBS')
    _Link('../../../../%s/frodan/%s/unbound/unbound.pdb' % (Env,strOfIndex), 'env.pdb')
    _Link('../../../../%s/frodan/%s/unbound/unbound.pdb' % (AB,strOfIndex), 'ab.pdb')
    _Link('../../../../../tools/VMD_mergePDB.tcl', 'VMD_mergePDB.tcl')
    cmdLine = 'vmd -e VMD_mergePDB.tcl -args env.pdb ab.pdb'
    rvals = subProcessCMD(cmdLine, 'Executing VMD_mergePDB.tcl for %s' % (Complex))
    if (rvals[2]):
        print 'PrepComplex - VMD returned error code: %s' % (rvals[2])
        writeLogs(rvals,'VMDmergePDB')
        raise Exception ('PrepComplex - VMD returned error code: %s' % (rvals[2]))
    myfiles = ['ab_fit.pdb', 'merged.ppsf', 'env_fit.pdb', 'seperated_ab.ali',\
                'seperated_ab.pap','seperated_env.ali','seperated_env.pap',\
                'seperated_fit.pdb']
    #for file in myfiles:
    #    _Remove(file)

def DockComplex(index,Complex):
    """Perturb the complex into the final docked state """
    structures = Complex.split('Complex-')[1]
    Env,AB = structures.split('___')
    strOfIndex = str(index).zfill(2)
    _SetWD('%s/Structures/%s/frodan/%s/bound' % (wDir,Complex,strOfIndex))
    if not(_Exists('bound.pdb')):
        _Link('../unbound/merged.pdb', 'merged.pdb')
        _Link('../../../../../tools/VMD_fixchainids.tcl', 'VMD_fixchainids.tcl')
        #_Link('../../../../../tools/VMD_align.tcl', 'VMD_align.tcl')
        _Link('../../../../../src/options.xml', 'options.xml')
        _Link('../../../../../src/grompp.mdp', 'grompp.mdp')
        _Link('../../../../../targets/%s' % (dockedTarget), 'docked.pdb')
        if(Global['multiChain'] == '1'):
            _Link('../../../../../src/Complexchains.txt', 'chains.txt')
            cmdLine = 'vmd -e VMD_fixchainids.tcl -args merged.pdb fixed.pdb 3 G H L'
        else:
            _Link('../../../../../src/ComplexCD4chains.txt', 'chains.txt')
            cmdLine = 'vmd -e VMD_fixchainids.tcl -args merged.pdb fixed.pdb 2 G C'
        rvals = subProcessCMD(cmdLine, 'Executing VMD_fixchainids.tcl for %s' % (Complex))
        if (rvals[2]):
            print 'DockComplex - VMD returned error code: %s' % (rvals[2])
            writeLogs(rvals,'VMDfixchainids')
            raise Exception('DockComplex - VMD returned error code: %s' % (rvals[2]))
        currentFile = 'fixed.pdb'
        _SetWD('%s/Structures/%s/frodan/%s/bound/' % (wDir,Complex,strOfIndex))
        if (multifrodan(currentFile, 'docked.pdb')):
            raise Exception('unrecoverable error in frodan',\
                    '%s/Structures/%s/frodan/%s/bound' % (wDir,Complex,strOfIndex))
        _Rename('frodan.xtc', 'bound.xtc')
        _Remove('conf.gro')
        prep_gmx('frodan01.pdb','%s/Structures/%s/frodan/%s/bound' % (wDir,Complex,strOfIndex))
        gmx_minimize('conf.gro','minimized_complex.pdb', \
                        '%s/Structures/%s/frodan/%s/bound/' % (wDir,Complex,strOfIndex))
        _SetWD('%s/Structures/%s/frodan/%s/bound/' % (wDir,Complex,strOfIndex))
        if(Global['multiChain'] == '1'):
            cmdLine = 'vmd -e VMD_fixchainids.tcl -args minimized_complex.pdb tmp.pdb 3 G H L'
        else:
            cmdLine = 'vmd -e VMD_fixchainids.tcl -args minimized_complex.pdb tmp.pdb 2 G C'
        rvals = subProcessCMD(cmdLine, 'Executing VMD_fixchainids.tcl for %s' % (Complex))
        if (rvals[2]):
            print 'DockComplex - VMD 2 returned error code: %s' % (rvals[2])
            writeLogs(rvals,'VMDfixchainids2')
            raise Exception('DockComplex - VMD 2 returned error code: %s' % (rvals[2]))
        try:
            salign('frodan01','tmp')
        except Exception as err:
            print 'DockComplex salign error - :', err
            raise Exception ('DockComplex salign error - :', err)
        _Rename('tmp_fit.pdb', 'bound.pdb')
        myfiles = ['conf.box.gro','confout.gro','fixed.pdb','tmp.pdb',\
                    'frodan00.pdb','frodan00.xtc','frodan01.pdb','#topol.tpr.1#',\
                    'minimized_model.pdb', 'model.pdb','md.log','mdout.mdp',\
                    'min_bound_complex.pdb', '#conf.box.gro.1#', '#confout.gro.1#',\
                    'posre_Protein_chain_G.itp','posre_Protein_chain_H.itp',\
                    'posre_Protein_chain_L.itp','topol_Protein_chain_G.itp',\
                    'topol_Protein_chain_H.itp','topol_Protein_chain_L.itp',\
                    'topol.top','topol.tpr','traj.trr','#conf.box.gro.1#',\
                    'minimized_complex.pdb','posre.itp','conf.gro','io.mc',\
                    'frodan01_fit.pdb','frodan01_tmp.ali','frodan01_tmp.pap',\
                    'distances_State1ToState2.txt']
        for thisFile in myfiles:
            if (_Exists(thisFile)):
                _Remove(thisFile)
        if(index%10 > 1):
            _Remove('bound.xtc')
# Undetermined use at this time
#TODO determine use of this method 
def SplitComplex(Complex,index,wDir):
    """Split the complex into bound Env and ABs and deposit into bound folders for each"""
    structures = Complex.split('Complex-')[1]
    Env,AB = structures.split('___')
    strOfIndex = str(index).zfill(2)
    if(Global['multiChain'] == '1'):
        stype = 'AB'
    else:
        stype = 'CD4'
    _SetWD('%s/Structures/%s/frodan/%s/bound' % (wDir,Complex,strOfIndex))
    _Link('../../../../../tools/VMD_sep_Complex.tcl', 'VMD_sep_Complex.tcl')
    cmdLine = ['vmd','-e','VMD_sep_Complex.tcl','-args', 'bound.pdb', stype,\
                '%s/Structures/%s/frodan/%s/bound/bound.pdb' % (wDir,Env,strOfIndex),\
                '%s/Structures/%s/frodan/%s/bound/bound.pdb' % (wDir,AB,strOfIndex)]
    ProcessCMD(cmdLine[0],cmdLine, 'Executing VMD_sep_Complex.tcl for %s' % (Complex))

def PrepFrodaN(model,seq,tag):
    """ working unit index and tag """
    frodanDir = '%s/Structures/%s/frodan/%s' % (wDir,seq,str(model).zfill(2))
    modellerDir = '%s/Structures/%s/modeller' % (wDir,seq)
    _SetWD(frodanDir)
    Write_grompp(frodanDir)
    if not(_Exists('model.pdb')):
        currentFile = '%s/ESSC.%s.pdb' % (modellerDir,str(model).zfill(4))
        prep_gmx(currentFile,frodanDir)
        gmx_minimize('conf.gro','model.pdb',frodanDir)
        _SetWD(frodanDir)
        if tag == 0:
            _Link('../../../../../targets/%s' % (boundTarget), '%s/bound/%s' % (frodanDir,boundTarget))
            _Link('../../../../../targets/%s' % (unboundTarget), '%s/unbound/%s' % (frodanDir,unboundTarget))
            filelist = ['posre.itp']
        elif tag == 1:
            _Rename('model.pdb','model.bak.pdb')
            _Link('../../../../tools/VMD_fixchainids.tcl','VMD_fixchainids.tcl')
            if(Global['multiChain'] == '1'):
                cmdLine = 'vmd -dispdev none -e VMD_fixchainids.tcl -args model.bak.pdb model.pdb 2 A B'
            else:
                cmdLine = 'vmd -dispdev none -e VMD_fixchainids.tcl -args model.bak.pdb model.pdb 1 A'
            rvals = subProcessCMD(cmdLine,'Executing VMD (chain label correction) for %s' % (currentFile))
            if (rvals[2]):
                print 'PrepFrodan - VMD returned error code: %s' % (rvals[2])
                writeLogs(rvals,'VMDfixchainids')
                raise Exception ('PrepFrodan - VMD returned error code: %s' % (rvals[2]))
            _Link('../../../../../targets/%s' % (boundABTarget), '%s/bound/%s' % (frodanDir,boundABTarget))
            _Link('../../../../../targets/%s' % (unboundABTarget), '%s/unbound/%s' % (frodanDir,unboundABTarget))
            if(Global['multiChain'] == '1'):
                filelist = ['topol_Protein_chain_A.itp','topol_Protein_chain_B.itp',\
                            'posre_Protein_chain_A.itp','posre_Protein_chain_B.itp']
            else:
                filelist = []
        filelist += ['topol.tpr','model.pdb','topol.top','conf.gro',\
                    'grompp.mdp','mdout.mdp','confout.gro']
        for state in ['bound','unbound']:
            _Link('../../../../../src/options.xml','%s/%s/options.xml' % (frodanDir,state))
            if tag == 0:
                _Link('../../../../../src/chains.txt','%s/%s/chains.txt' % (frodanDir,state))
            elif tag == 1:
                if(Global['multiChain'] == '1'):
                    _Link('../../../../../src/ABchains.txt','%s/%s/chains.txt' % (frodanDir,state))
                else:
                    _Link('../../../../../src/CD4chains.txt','%s/%s/chains.txt' % (frodanDir,state))
            for fileName in filelist:
                _Link('../%s' % (fileName),'%s/%s/%s' % (frodanDir,state,fileName))

def ProcessCMD(cmd,args, message):
    if len(message) > 0:
        print(message)
    thisDir = os.getcwd()   
    rc = os.fork()
    if rc == -1:
        print 'fork failed'
        raise OSError('os.fork() failed')
    elif rc == 0:
        _SetWD(thisDir)
        sys.stdout.flush()
        try:
            os.execvp(cmd,args)
        except Exception as err:
            print '---Child process unable to execute command "%s"\n, Error is-%s' % (cmd,err)
    else:
        status = os.waitpid(rc,0)
        if not(status[1] == 0):
            sys.exit(-1)
            #return status[1]
        else:
            print 'Command %s with message %s completed\n' % (cmd,message)
            return 0

def subProcessCMD(cmd, message):
    rvals = []
    if len(message) > 0:
        print(message)
    try:
        p = sp.Popen(cmd,stdout=sp.PIPE,stderr=sp.PIPE,shell=True)
        rvals = list(p.communicate())
        rvals.append(p.wait())
    except Exception as err:
        print 'subProcessCMD error: {0}'.format(err)
        sys.exit(-1)
    print 'Command: "%s" completed\n' % (cmd)
    return rvals

def procTMPDX(args):
    targetDir,DX = args
    Dx = _ReadFile(DX)
    try:
        thisData = Dx[Dx.find('\n',Dx.rfind('data follows\n'))+1:Dx.find('attribute')].split()
    except Exception as  err:
        print err
        thisData = []
    if Global['CompressData'] == '1':
        if(DX.find('Complex-') > 0):
            tmpFile = DX.split('%s/' % (Global['tempDir']))[1].replace('_essc_','_')
            DXout = '%s/%s.ZDX' % (targetDir,tmpFile)
            DXmeta = '%s/%s' % (targetDir,tmpFile)
        else:
            DXout = '%s/%s.dx.ZDX' % (targetDir,targetDir[targetDir.rfind('/')+1:])
            DXmeta = '%s/%s.dx' % (targetDir,targetDir[targetDir.rfind('/')+1:])
        print 'Must feed the Cython monster hahahahaha'
        binArray = np.asarray(thisData,dtype=np.float64)
        header_footer = Dx[:Dx.find('\n',Dx.rfind('data follows\n'))]
        header_footer += '\n\n---\n%s.dx.ZDX\n---\n\n' % (targetDir[targetDir.rfind('/')+1:])
        header_footer += Dx[Dx.find('attribute'):]
        z,y,x = Dx[Dx.find('ts ',Dx.find('gridpositions counts'))+3:Dx.find('\n',Dx.find('gridpositions counts'))].split()
        args = '-b 1 -z %s -d -3 %s %s %s -a %s' % (DXout,x,y,z,Global['ZFPMaxError'])
        try:
            print 'Feeding the Cython monster... so hungry it is'
            rc = zfp.ZFPC(args,binArray)
        except Exception as err:
            print 'ZFP error:', err
            sys.exit(-1)
        if not(rc):
            try:
                thisFile = open(DXmeta,'w')
                thisFile.write(header_footer)
                thisFile.flush()
                os.fsync(thisFile)
                thisFile.close()
                _Remove(DX)
            except Exception as  err:
                _Remove(DX)
                print 'Post ZFP error:', err
                sys.exit(-1)
        else:
            print 'ZFP error\n'
    else:
        DXout = '%s/%s.dx' % (targetDir,targetDir[targetDir.rfind('/'):])
        _cp(DX,DXout)
        _Remove(DX)
    return thisData

def ReadpHFile(absfileName):
    """ Used to read and return pH/pH-conc-grid-subsets in list format"""
    try:
        thisFile = open(absfileName,'r')
        rtrnList = thisFile.read().splitlines()
        thisFile.close()
    except Exception as  err:
        print err
    return rtrnList
def reCalcPOTs(argsList):
    """ Single execution block"""
    targetDir,seq, state, model, pH = argsList
    rows = 0
    cols = 0    
    PQRSIZE = '%s/%s.size' % (targetDir,targetDir[targetDir.rfind('/'):])
    oldPotFile = np.loadtxt('%s.pot' % (PQRSIZE[:PQRSIZE.rfind('.')]))
    if not(np.max(oldPotFile) == 0.0 and np.min(oldPotFile) == 0.0):
        return
    DX = '%s/%s.dx' % (targetDir,targetDir[targetDir.rfind('/')+1:])
    SAS = '%s.sas.bin' % (PQRSIZE[:PQRSIZE.rfind('.')])
    Size = _ReadFile(PQRSIZE)
    inVars = ['FINEGRID','FINEDIM','GRIDCENTER']
    ofInterest = ['Num. fine grid pts. = ','Fine grid dims = ','Center = ']
    for i in range(len(inVars)):
        indexL = Size.find(ofInterest[i])
        indexU = Size.find('\n',indexL)
        inVars[i] = Size[indexL + len(ofInterest[i]):indexU].replace('x ','').replace('A','').rstrip().split()
    temp = np.linspace(0,float(inVars[1][0]),float(inVars[0][0]))
    temp = temp - np.mean(temp)
    tempX = temp + float(inVars[2][0])
    temp = np.linspace(0,float(inVars[1][1]),float(inVars[0][1]))
    temp = temp - np.mean(temp)
    tempY = temp + float(inVars[2][1])
    temp = np.linspace(0,float(inVars[1][2]),float(inVars[0][2]))
    temp = temp - np.mean(temp)
    tempZ = temp + float(inVars[2][2])
    mygrid = [list(tempX),list(tempY),list(tempZ)]
    mydiff = [tempX[1]-tempX[0],tempY[1]-tempY[0],tempZ[1]-tempZ[0]]
    try:
        Sas = np.fromfile(SAS)
        rows = np.shape(Sas)[0]
        Sas = np.reshape(Sas,(rows/3,3))
        rows,cols = np.shape(Sas)
        mysasInd = np.zeros((rows,cols))
    except Exception as  err:
        print ('reCalcSASPots - load SAS: {0}'.format(err))
    for myk in range(3):
        mysasSorted = np.sort(Sas[:,myk].view())
        mysasSortedix = Sas[:,myk].argsort(0)
        mysasPos = np.zeros(rows)
        myi = 0
        for myj in range(rows):
            while (myi < len(mygrid) and mysasSorted[myj] >= mygrid[myk][myi]):
                myi = myi + 1
            mysasPos[myj] = (myi + 1) - ((mygrid[myk][myi] - mysasSorted[myj])/mydiff[myk])
        mysasInd[mysasSortedix,myk] = mysasPos
    def mytranslate(x,y):
        index = (((x[:,0] - 1) * len(y[1]) * len(y[2])) + ((x[:,1] - 1) * len(y[2])) + x[:,2]) - 1
        return index.astype(int)
    mysasLow = mysasInd - np.floor(mysasInd)
    mysasHigh = np.ceil(mysasInd) - mysasInd
    results = np.zeros((rows,cols))
    Dx = _ReadFile(DX)
    if Global['CompressData'] == '1':
        DXin = '%s/%s.dx.ZDX' % (targetDir,targetDir[targetDir.rfind('/')+1:])
        z,y,x = Dx[Dx.find('ts ',Dx.find('gridpositions counts'))+3:Dx.find('\n',Dx.find('gridpositions counts'))].split()
        mydata = np.zeros((int(x),int(y),int(z)),dtype=np.float64) # I must feed the Cython monster
        args = '-b 0 -z %s -d -3 %s %s %s -a %s' % (DXin,x,y,z,Global['ZFPMaxError'])
        mydata = zfp.ZFPD(args,mydata)
        mydata = np.reshape(mydata,(int(x)*int(y)*int(z)))
    else:
        pass
        #TODO code up for an raw dx file
    try:
        TASKS = [np.asfarray(np.take(mydata,mytranslate(np.column_stack((np.floor(mysasInd[:,0]),\
                                                                                np.floor(mysasInd[:,1]),\
                                                                                np.floor(mysasInd[:,2])))\
                                                                                ,mygrid)))[:,None] * \
                                                np.column_stack((mysasLow[:,0],mysasLow[:,1],mysasLow[:,2])),\
                np.asfarray(np.take(mydata,mytranslate(np.column_stack((np.ceil(mysasInd[:,0]),\
                                                                                np.floor(mysasInd[:,1]),\
                                                                                np.floor(mysasInd[:,2])))\
                                                                                ,mygrid)))[:,None] * \
                                                np.column_stack((mysasHigh[:,0],mysasLow[:,1],mysasLow[:,2])),\
                                                
                np.asfarray(np.take(mydata,mytranslate(np.column_stack((np.ceil(mysasInd[:,0]),\
                                                                                np.ceil(mysasInd[:,1]),\
                                                                                np.floor(mysasInd[:,2])))\
                                                                                ,mygrid)))[:,None] * \
                                                np.column_stack((mysasHigh[:,0],mysasHigh[:,1],mysasLow[:,2])),\
                np.asfarray(np.take(mydata,mytranslate(np.column_stack((np.ceil(mysasInd[:,0]),\
                                                                                np.ceil(mysasInd[:,1]),\
                                                                                np.ceil(mysasInd[:,2])))\
                                                                                ,mygrid)))[:,None] * \
                                                np.column_stack((mysasHigh[:,0],mysasHigh[:,1],mysasHigh[:,2])),\
                np.asfarray(np.take(mydata,mytranslate(np.column_stack((np.floor(mysasInd[:,0]),\
                                                                                np.ceil(mysasInd[:,1]),\
                                                                                np.ceil(mysasInd[:,2])))\
                                                                                ,mygrid)))[:,None] * \
                                                np.column_stack((mysasLow[:,0],mysasHigh[:,1],mysasHigh[:,2])),\
                np.asfarray(np.take(mydata,mytranslate(np.column_stack((np.floor(mysasInd[:,0]),\
                                                                                np.floor(mysasInd[:,1]),\
                                                                                np.ceil(mysasInd[:,2])))\
                                                                                ,mygrid)))[:,None] * \
                                                np.column_stack((mysasLow[:,0],mysasLow[:,1],mysasHigh[:,2])),\
                np.asfarray(np.take(mydata,mytranslate(np.column_stack((np.floor(mysasInd[:,0]),\
                                                                                np.ceil(mysasInd[:,1]),\
                                                                                np.floor(mysasInd[:,2])))\
                                                                                ,mygrid)))[:,None] * \
                                                np.column_stack((mysasLow[:,0],mysasHigh[:,1],mysasLow[:,2])),\
                np.asfarray(np.take(mydata,mytranslate(np.column_stack((np.ceil(mysasInd[:,0]),\
                                                                                np.floor(mysasInd[:,1]),\
                                                                                np.ceil(mysasInd[:,2])))\
                                                                                ,mygrid)))[:,None] * \
                                                np.column_stack((mysasHigh[:,0],mysasLow[:,1],mysasHigh[:,2]))]
        pool = Pool(4)        
        resultList = pool.imap(calculate,TASKS)
        pool.close()
        pool.join()
        for r in resultList:
            results = np.add(results, r)
    except Exception as  err:
        print argsList
        print 'reCalcSASPots (average charge at a point) : {0}'.format(err)
        sys.exit(-1)
    results = np.mean(results / 4,axis=1)
    try:    
        np.savetxt('%s.pot' % (PQRSIZE[:PQRSIZE.rfind('.')]),results,fmt='%12.10f')
    except Exception as  err:
        print 'recalcPots error :', err

def RunFrodaN(index,seq,state,tag):
    """ as is """
    frodanDir = '%s/Structures/%s/frodan/%s/%s' % (wDir,seq,str(index).zfill(2),state)
    _SetWD(frodanDir)
    _Link('../../../../../tools/VMD_fixchainids.tcl', 'VMD_fixchainids.tcl')
    if tag == 0:
        bTarget = boundTarget
        ubTarget = unboundTarget
    elif tag == 1:
        bTarget = boundABTarget
        ubTarget = unboundABTarget
    if (state =='bound'):
        if not(_Exists('bound.pdb')):
            if (multifrodan('model.pdb', bTarget)):
                print 'unrecoverable error in frodan, creating new model'
                recList = copy.copy(argsList)
                Get_New_Model(recList,tag)
                _SetWD(frodanDir)
                if (multifrodan('model.pdb', bTarget)):
                    raise Exception('unrecoverable error in frodan',frodanDir)
            _Rename('frodan.xtc', 'bound.xtc')
            _Remove('conf.gro')
            gmx_minimize('frodan01.pdb','bound.pdb',frodanDir)
            _SetWD(frodanDir)
            if(Global['antibody'] == '1'):
                if(abseqs.__contains__(seq)):
                    _Rename('bound.pdb','ab_tmp.pdb')
                    _Link('../../../../../targets/%s' % (ABsalign), 'align_sep.pdb')
                    if(Global['multiChain'] == '1'):
                        cmdLine = 'vmd -e VMD_fixchainids.tcl -args ab_tmp.pdb tmp.pdb 2 H L'
                    else:
                        cmdLine = 'vmd -e VMD_fixchainids.tcl -args ab_tmp.pdb tmp.pdb 1 C'
                    rvals = subProcessCMD(cmdLine, 'Executing VMD_fixchainids.tcl for %s' % (seq))
                    if (rvals[2]):
                        writeLogs(rvals,'VMDfixchainids')
                        raise Exception ('RunFrodaN - VMD (bound) returned error code: %s' % (rvals[2]))
                else:
                    _Rename('bound.pdb','tmp.pdb')
                    _Link('../../../../../targets/%s' % (gp120salign), 'align_sep.pdb')
                try:
                    salign('align_sep','tmp')
                except Exception as err:
                    print 'RunFrodan - salign error (bound) :', err
                _Rename('tmp_fit.pdb', 'bound.pdb')
    elif (state =='unbound'):
        if not(_Exists('unbound.pdb')):
            if(multifrodan('model.pdb', ubTarget)):
                print 'unrecoverable error in frodan, creating new model'
                recList = copy.copy(argsList)
                Get_New_Model(recList,tag)
                _SetWD(frodanDir)
                if(multifrodan('model.pdb', ubTarget)):
                    raise Exception('unrecoverable error in frodan',frodanDir)
            _Rename('frodan.xtc', 'unbound.xtc')
            _Remove('conf.gro')
            gmx_minimize('frodan01.pdb','unbound.pdb',frodanDir)
            _SetWD(frodanDir)
            if(Global['antibody'] == '1'):
                if(abseqs.__contains__(seq)):
                    _Rename('unbound.pdb','ab_tmp.pdb')
                    _Link('../../../../../targets/%s' % (ABsalign), 'align_sep.pdb')
                    if(Global['multiChain'] == '1'):
                        cmdLine = 'vmd -e VMD_fixchainids.tcl -args ab_tmp.pdb tmp.pdb 2 H L'
                    else:
                        cmdLine = 'vmd -e VMD_fixchainids.tcl -args ab_tmp.pdb tmp.pdb 1 C'
                    rvals = subProcessCMD(cmdLine, 'Executing VMD_fixchainids.tcl for %s' % (seq))
                    if (rvals[2]):
                        writeLogs(rvals,'VMDfixchainids')
                        raise Exception ('RunFrodaN - VMD (unbound) returned error code: %s' % (rvals[2]))
                else:
                    _Rename('unbound.pdb','tmp.pdb')
                    _Link('../../../../../targets/%s' % (gp120salign), 'align_sep.pdb')
                try:
                    salign('align_sep','tmp')
                except Exception as err:
                    print 'RunFrodan - salign error (unbound) :', err
                _Rename('tmp_fit.pdb', 'unbound.pdb')
    else:
        print('Error: unknown state provided: ',state)
    for file in glob('frodan*'):
        _Remove(file)
    myfiles = ['#confout.gro.1#','#mdout.mdp.1#','conf.box.gro','confout.gro','md.log',\
                'ener.edr','traj.trr','#topol.tpr.1#','tmp.pdb','align_sep_fit.pdb',\
                'mdout.mdp','topol.tpr','align_sep_tmp.ali','align_sep_tmp.pap']
    for thisFile in myfiles:
        if (_Exists(thisFile)):
            _Remove(thisFile)
    if(index%10 > 1):
        _Remove('%s.xtc' % (state))

def Run_Modeller(seq,tag):
    """ Runs modeller against the specified sequence file """
    modellerDir = '%s/Structures/%s/modeller' % (wDir,seq)
    startingModel = 1
    _SetWD(modellerDir)
    pdbs = glob('ESSC.*.pdb')
    if(len(pdbs) == numModels):
        return
    else:
        startingModel = numModels - len(pdbs)
    if tag == 0:
        for fileName in glob('../../../templates/template?.pdb'):
            _Link(fileName,fileName.split('/')[-1])
        Convert_Seq('../../../seqfiles/%s.seq' % (seq), 'sequence.ali')
    else:
        for fileName in glob('../../../ABtemplates/template?.pdb'):
            _Link(fileName,fileName.split('/')[-1])
        Convert_Seq('../../../ABseqfiles/%s.seq' % (seq), 'sequence.ali')
    print 'Modeller is running at: %s' % (modellerDir)
    #homomodel(numModels)
    hmodel = homomodel.replace('++!!SART_MODEL_NUM!!++',str(startingModel))
    hmodel = homomodel.replace('++!!NUM_MODELS!!++',str(numModels))
    thisFile = open('homomodel.py','w')
    thisFile.write(hmodel)
    thisFile.flush()
    os.fsync(thisFile)
    thisFile.close()
    cmdLine = 'chmod +x homomodel.py'
    rvals = subProcessCMD(cmdLine,'Setting execution bit for homomodel.py')
    if (rvals[2]):
        writeLogs(rvals,'chmod')
        raise Exception ('Modeller - chmod died setting the execution bit, must have been wearing a red shirt')
    cmdLine = 'modpy.sh ./homomodel.py'
    rvals = subProcessCMD(cmdLine,'Executing Modeller')
    writeLogs(rvals,'modeller')

def salign(target,fitting):
    print 'salign for target : fitting ', target, fitting
    #log.minimal()
    log.level(output=0, notes=0, warnings=0, errors=1, memory=0)
    log.none()
    env = environ()
    env.io.atom_files_directory = ['./','../atom_files/']
    aln = alignment(env)
    codes=[]
    ## If you have multiple chains, multiple templates...
    for (code) in (target,fitting):
        codes.append(code)
        mdl = model(env, file='%s.pdb' % (code))
        aln.append_model(mdl, atom_files=code, align_codes=code)
    fName = '%s_%s' % (codes[0],codes[1])
    for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True),
                                        ((1., 0.5, 1., 1., 1., 0.), False, True),
                                        ((1., 1., 1., 1., 1., 0.), True, False)):
        aln.salign(rms_cutoff=3.5, normalize_pp_scores=False,
                   rr_file='$(LIB)/as1.sim.mat', overhang=30,
                   gap_penalties_1d=(-450, -50),
                   gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0,
                   dendrogram_file='%s_fm00495.tree' % (fName),
                   alignment_type='tree', # If 'progresive', the tree is not
                                          # computed and all structues will be
                                          # aligned sequentially to the first
                   feature_weights=weights, # For a multiple sequence alignment only
                                            # the first feature needs to be non-zero
                   improve_alignment=True, fit=True, write_fit=write_fit,
                   write_whole_pdb=whole, output='ALIGNMENT QUALITY',
                   edit_file_ext=('.pdb','_fit.pdb'))
    aln.write(file='%s.pap' % (fName), alignment_format='PAP')
    aln.write(file='%s.ali' % (fName), alignment_format='PIR')
    aln.salign(rms_cutoff=1.0, normalize_pp_scores=False,
               rr_file='$(LIB)/as1.sim.mat', overhang=30,
               gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3),
               gap_gap_score=0, gap_residue_score=0, 
               dendrogram_file='%s_1is3A.tree' % (fName),
               alignment_type='progressive', feature_weights=[0]*6,
               improve_alignment=False, fit=False, write_fit=True,
               write_whole_pdb=False, output='QUALITY')
    print 'Finished salign for target : fitting ', target, fitting

def Validate_CFG_File(config):
    #Check the config file for errors or missing parameters
    lexica = GetLanguage()
    errorFlag = False
    errors = []
    ZFPError = lexica['Global']['ZFPError']
    studies = config['Studies']
    globalE = config['Global']
    for name in globalE:
        if name == 'tempDir':
            try:
                os.mkdir('%s/essc' % (globalE['tempDir']))
                os.rmdir('%s/essc' % (globalE['tempDir']))
            except Exception as  err:
                errorFlag = True
                errors.append("Error in Global: '%s' unable to write to temp directory %s" % (name,str(err)))
        elif name == 'ZFPMaxError':
            if not(ZFPError.__contains__(globalE[name])):
                errorFlag = True
                errors.append("Error in Global: '%s' does not exist or is not a valid entry" % (name))
        elif (not(lexica['Global']['values'].__contains__(name)) or not(globalE[name].isdigit())):
            errorFlag = True
            errors.append("Error in Global: entry '%s' does not exist or is not a number" % (name))
    dirs = lexica['Studies']['dirs']
    files = lexica['Studies']['files']
    values = lexica['Studies']['values']
    for key in studies.keys():
        current = studies[key]
        for name in dirs:
            if (not(current.has_key(name)) or not(_Exists(studies[key][name]))):
                errorFlag = True
                errors.append("Error in seq: {%s} entry or directory for '%s' does not exist" % (key,name))
        for name in files:
            if (not(current.has_key(name))):
                errorFlag = True
                errors.append("Error in seq: {%s} entry does not exist" % (key))
            if (name == 'boundTarget' or name == 'unboundTarget' or name == 'dockedTarget'):
                    if not(_Exists('%s/targets/%s' % (studies[key]['wdir'],studies[key][name]))):
                        errorFlag = True
                        errors.append("Error in seq: {%s} file for '%s' does not exist" % (key,name))
            if (name == 'cgridsubsetFile'):
                    if not(_Exists(studies[key][name])):
                        errorFlag = True
                        errors.append("Error in seq: {%s} file for '%s' does not exist" % (key,name))
        for name in values:
            if (not(current.has_key(name)) or not(studies[key][name].isdigit())):
                errorFlag = True
                errors.append("Error in seq: {%s} entry '%s' does not exist or is not a number" % (key,name))
        if not(_Exists('%s/src/options.xml' % (studies[key]['wdir']))):
            errorFlag = True
            errors.append('Error: src/options.xml file does not exist')
        if not(_Exists('%s/src/chains.txt' % (studies[key]['wdir']))):
            errorFlag = True
            errors.append('Error: src/chains.txt file does not exist')
        if not(_Exists('%s/src/homomodel.py' % (studies[key]['wdir']))):
            errorFlag = True
            errors.append('Error: src/homomodel.py does not exist')
        if not(_Exists('%s/src/generic.in' % (studies[key]['wdir']))):
            errorFlag = True
            errors.append('Error: src/generic.in does not exist')
        if not(globalE['SMT-HT'] == '1' or globalE['SMT-HT'] == '2'):
            errorFlag = True
            errors.append('Error: Global:SMT-HT not valid value')
        if (globalE['ABmode'] =='1'):
            if not(_Exists('%s/tools/VMD_mergePDB.tcl' % (studies[key]['wdir']))):
                errorFlag = True
                errors.append('Error: tools/VMD_mergePDB.tcl does not exist')
            if not(_Exists('%s/tools/VMD_merge_all_Complex.tcl' % (studies[key]['wdir']))):
                errorFlag = True
                errors.append('Error: tools/VMD_merge_all_Complex.tcl does not exist')
            if not(_Exists('%s/tools/VMD_fixchainids.tcl' % (studies[key]['wdir']))):
                errorFlag = True
                errors.append('Error: tools/VMD_fixchainids.tcl does not exist')
            if not(_Exists('%s/src/grompp.mdp' % (studies[key]['wdir']))):
                errorFlag = True
                errors.append('Error: src/grompp.mdp does not exist')
            if not(_Exists('%s/src/skipAPBS' % (studies[key]['wdir']))):
                errorFlag = True
                errors.append('Error: src/skipAPBS does not exist')
            if not(_Exists('%s/src/complex.in' % (studies[key]['wdir']))):
                errorFlag = True
                errors.append('Error: src/complex.in does not exist')
            if(globalE['multiChain'] == '1'):
                if not(_Exists('%s/src/ABchains.txt' % (studies[key]['wdir']))):
                    errorFlag = True
                    errors.append('Error: src/ABchains.txt file does not exist')
                if not(_Exists('%s/src/Complexchains.txt' % (studies[key]['wdir']))):
                    errorFlag = True
                    errors.append('Error: src/Complexchains.txt file does not exist')
            else:
                if not(_Exists('%s/src/CD4chains.txt' % (studies[key]['wdir']))):
                    errorFlag = True
                    errors.append('Error: src/CD4chains.txt file does not exist')
                if not(_Exists('%s/src/ComplexCD4chains.txt' % (studies[key]['wdir']))):
                    errorFlag = True
                    errors.append('Error: src/ComplexCD4chains.txt file does not exist')
        if(int(globalE['doRAnalysis'])):
            if not(_Exists('%s/BESI_EVM_analysis.r' % (studies[key]['wdir']))):
                errorFlag = True
                errors.append('Error: BESI_EVM_analysis.r does not exist')
            if not(_Exists('%s/libBESI_EVM.r' % (studies[key]['wdir']))):
                errorFlag = True
                errors.append('Error: libBESI_EVM.r does not exist')
    if errorFlag:
        for entry in errors:
            print (entry)
        sys.exit(-1)
    else:
        print('No errors in configuration or required files\nProceeding with run')
def Write_grompp(dest):
    _Link('../../../../src/grompp.mdp','%s/grompp.mdp' % (dest))
def writeLogs(data,fileName):
    oFile = gzip.open('%s.log.gz' % fileName,'wb')
    oFile.write(data[0])
    oFile.close()
    if (len(data[1]) > 0):
        oFile =  gzip.open('%s.err.gz' % fileName,'wb')
        oFile.write(data[1])
        oFile.close()
# Init the lib
__init__()
