#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    ESSC pipeline driver <pline.py>
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
#   Name:pline_Driver.py  
#   Python Script
#   Author: Scott P Morton (spm3c at mtmail.mtsu.edu)
# 
#   pline_Driver is an MPI based script for the essc-pipeline library
#   This version is specifically developed for Anti-Body studies
#   where bnAbs are the focus for Complex targets. For instance,
#   b12 is bound to 2NY7 at CD4. We will use 2NY7 as the target for 
#   the complex conformation
#
#
#   Requires 'libessc.py'
#
#---------------------------
"""

__title__  = 'pline'
__author__ = 'Scott P Morton (spm3c at mtmail.mtsu.edu)'
Description = '\
    Driver for the Electro-Static Surface Charge Pipeline Library\n\n \
        Copyright (C) 2018  MIDDLE TENNESSEE STATE UNIVERSITY\n \
        Copyright (C) 2018  Scott P. Morton\n \
        \n\n \
         __  __ _____ ____  _   _ \n \
        |  \/  |_   _/ ___|| | | |\n \
        | |\/| | | | \___ \| | | |\n \
        | |  | | | |  ___) | |_| |\n \
        |_|  |_| |_| |____/ \___/ \n \
        \n\n \
        MIDDLE TENNESSEE STATE UNIVERSITY\n\n \
        This program comes with ABSOLUTELY NO WARRANTY;\n \
        This is free software, and you are welcome to redistribute it\n \
        under certain conditions; see <http://www.gnu.org/licenses/> for details.\n\n'

def Version():
    return '1.1.6 pline.py'

import sys, copy, time, os
import libessc as essc

#import testlib as essc

import numpy as np
from mpi4py import rc
rc.initialize = False
  
from mpi4py import MPI
assert not MPI.Is_initialized()
assert not MPI.Is_finalized()
  
MPI.Init_thread(MPI.THREAD_MULTIPLE)
assert MPI.Is_initialized()
assert not MPI.Is_finalized()

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
world_Size = comm.Get_size()
thisHost = os.getenv('HOSTNAME')
if rank == 0:
    print Version()
    print Description
    essc.printInfo()

# Variables
seqFiles = []
ABseqFiles = []
myseqFiles = []
myABseqFiles = []
allSeqFiles = []
complexStructures = []
mycomplexStructures = []
allStructures = []
meansIndexList = []
cgridsubset = []
homomodel = ''
complexIN = ''
genericIN = ''

config = {}
myconfig = {}

if rank == 0:
    print 'Setting up environment'
    print time.ctime()
    if(len(sys.argv) > 1):
        if(sys.argv[1].upper() == '--HELP' or sys.argv[1].upper() == '-H'):
            with open('README','r') as f:
                readme = f.read()
            print readme
        else:
            config = essc._ReadConfig(sys.argv[1])
            wDir = config['Studies']['bulk_mode']['wdir']
            gridFile = config['Studies']['bulk_mode']['cgridsubsetFile']
            config['Studies']['bulk_mode']['templates'] = '%s/templates' % (wDir)
            config['Studies']['bulk_mode']['seqFiles'] = '%s/seqfiles' % (wDir)
            config['Studies']['bulk_mode']['sourceFiles'] = '%s/src' %(wDir)
            config['Studies']['bulk_mode']['cgridsubsetFile'] = '%s/src/%s' % (wDir,gridFile)
            config['Studies']['bulk_mode']['ABseqFiles'] = '%s/ABseqfiles' % (wDir)
            config['Studies']['bulk_mode']['ABtemplates'] = '%s/ABtemplates' % (wDir)
            #Check the environment for errors or missing parameters/files
            essc.Validate_CFG_File(config)
            essc.Global=config['Global']
            essc.Studies=config['Studies']
            Files = essc._ReadSeqFiles(config['Studies']['bulk_mode']['seqFiles'])
            for i in range(len(Files)):
                seqFiles.append(Files[i][Files[0].rfind('/')+1:Files[i].rfind('.')])
            if config['Global']['antibody'] == '1':
                # Driver must load and broadcast AB files and APBS .in files
                ABFiles = essc._ReadSeqFiles(config['Studies']['bulk_mode']['ABseqFiles'])
                genericIN = essc._ReadFile('%s/src/generic.in' % (wDir))
                complexIN = essc._ReadFile('%s/src/complex.in' % (wDir))
                for i in range(len(ABFiles)):
                    AB = ABFiles[i][ABFiles[0].rfind('/')+1:ABFiles[i].rfind('.')]
                    ABseqFiles.append(AB)
                if config['Global']['ABmode'] == '0':
                    for env in seqFiles:
                        for ab in ABseqFiles:
                            complexStructures.append('Complex-%s___%s' % (env,ab))
                elif config['Global']['ABmode'] == '1':
                    abmaps = config['ABmaps']
                    for key in abmaps.keys():
                        complexStructures.append('Complex-%s___%s' % (key,abmaps[key]))
            cgridsubset = essc.ReadpHFile('%s/src/%s' % (wDir,gridFile))
            homomodel = essc._ReadFile('%s/src/homomodel.py' % (wDir))

comm.Barrier()

# Broadcast base requirements to all ranks
myconfig = comm.bcast(config)
myseqFiles = comm.bcast(seqFiles)
essc.cgridsubsetList = comm.bcast(cgridsubset)
essc.homomodel = comm.bcast(homomodel)
essc.genericIN = comm.bcast(genericIN)

# Setup AB specific needs and broadcast
if myconfig['Global']['antibody'] == '1':
    myABseqFiles = comm.bcast(ABseqFiles)
    mycomplexStructures = comm.bcast(complexStructures)
    essc.complexIN = comm.bcast(complexIN)

allSeqFiles = myseqFiles + myABseqFiles
allStructures = allSeqFiles + mycomplexStructures

# Set library global variables !!!! MUST HAVE !!!!
essc.abseqs = myABseqFiles
essc.Global=myconfig['Global']
essc.Studies=myconfig['Studies']
essc.myrank = rank
#seq_env=myconfig['Studies']
if myconfig.has_key('ABmaps'):
    essc.ABmaps=myconfig['ABmaps']

# Finish initing the lib
essc.completeInit()

# Setup global variables for repeated use
numModels = essc.numModels
limitAPBS = essc.limitAPBS
states = essc.states
pHgridLen = essc.pHgridLen
lenCmplxStrs = len(mycomplexStructures)
lenAllSeqFiles = len(allSeqFiles)

# determine limiting value for certain operations
if (numModels > limitAPBS):
    numMods = limitAPBS
else:
    numMods = numModels


MPIbuffer=np.zeros(1,dtype=np.int64)
workIDX = np.int64(0)
ESSCstatus = MPI.Status()

essc._SetWD(essc.wDir)

if rank == 0:
    # Rank 0 is the work scheduler
    messages = ['Preparing to build directories for Complex Structures',\
                'Preparing to build directory structures and run Modeller',\
                'Setting up to prep FrodaN',\
                'Setting up to run FrodaN',\
                'Setting up to run PrepComplex',\
                'Setting up to run DockComplex',\
                'Setting up to run MergeALLPDB',\
                'Setting up to run GenPsize',\
                'Setting up to run genComplexIN',\
                'Setting up to run genIN',\
                'Setting up to run GenSAS processes',\
                'Setting up to run APBS processes',\
                'Setting up to calculate means and ratios for all individual structures',\
                'Setting up to calculate Binding Energies for Complex structures',\
                'Setting up to calculate Residue surface charges',\
                'Setting up to correlate Residue surface charge data']
    # Generate the work units for each section of work
    genpsizeWork = np.arange(0,1,dtype=np.int64)
    structuresWork = np.arange(0,lenCmplxStrs,1,np.int64)
    modellerWork = np.arange(0,lenAllSeqFiles,1,np.int64)
    prepFrodaNWork = np.arange(0,lenAllSeqFiles * numModels,1,np.int64)
    runFrodaNWork = np.arange(0,lenAllSeqFiles * numModels * 2,1,np.int64)
    ComplexWork = np.arange(0,lenCmplxStrs * numModels,1,np.int64)
    ComplexgenINWork = np.arange(0,lenCmplxStrs * numModels * pHgridLen,1,np.int64)
    APBSWork = np.arange(0,len(allStructures) * numModels * pHgridLen * 2,1,np.int64)
    APBSWorkLTD = np.arange(0,len(allStructures) * numMods * pHgridLen * 2,1,np.int64)
    ResWork = np.arange(0,lenAllSeqFiles * numMods * pHgridLen * 2,1,np.int64)
    # ResWork = np.arange(0,len(myseqFiles) * numMods * pHgridLen * 2,1,np.int64)
    resultsWork = np.arange(0,lenAllSeqFiles,1,np.int64)
    resMnRWork = np.arange(0,lenAllSeqFiles * 2,1,np.int64)
    Stages = [structuresWork,modellerWork,
              prepFrodaNWork,runFrodaNWork,
              ComplexWork,ComplexWork,
              ComplexWork,APBSWork,
              ComplexgenINWork,APBSWork,APBSWork,APBSWorkLTD,
              resultsWork,structuresWork,
              ResWork,resMnRWork]
    workType =  range(0,len(Stages))
    #workType =  range(13,14) # For stepping through
    for wtype in workType:
        print messages[wtype]
        print time.ctime()
        for idx in Stages[wtype]:
            comm.Recv(MPIbuffer, source=MPI.ANY_SOURCE,tag=10,status=ESSCstatus)
            srcRank = ESSCstatus.Get_source()
            ESSCrequest = comm.Isend(idx,dest=srcRank,tag=wtype)
            ESSCrequest.Wait(ESSCstatus)
        for i in np.arange(1,world_Size,1,np.int64):
            comm.Recv(MPIbuffer, source=MPI.ANY_SOURCE,tag=10,status=ESSCstatus)
            srcRank = ESSCstatus.Get_source()
            ESSCrequest = comm.Isend(i,dest=srcRank,tag=100)
            ESSCrequest.Wait(ESSCstatus)
        comm.Barrier()
    for i in np.arange(1,world_Size,1,np.int64):
        comm.Recv(MPIbuffer, source=MPI.ANY_SOURCE,tag=10,status=ESSCstatus)
        srcRank = ESSCstatus.Get_source()
        ESSCrequest = comm.Isend(i,dest=srcRank,tag=200)
        ESSCrequest.Wait(ESSCstatus)
else:
    # Build basic reuse variables for task operations
    tag=0
    messages = ['Build Directories for Complex structures',\
                'Build directory structures and run Modeller',\
                'prep FrodaN',\
                'Run FrodaN',\
                'PrepComplex',\
                'DockComplex',\
                'Merge PDBs for Complex',\
                'GenPsize for All Structures',\
                'Run genComplexIN',\
                'Run genIN',\
                'Run GenSAS',\
                'Run APBS',\
                'Calculate means and ratios for base structures',\
                'Calculate Complex Structure Binding Energies',\
                'Calculating residue surface charges',\
                'Correlating residue surface charges']
    #### !!!!!!!!!!!!!!! ####
    # Define Pipeline
    def Complex():
        seq = copy.copy(mycomplexStructures[workIDX])
        essc.CreateDirectoryStructure(seq)
    def modeller():
        myTag = 0
        if workIDX >= len(myseqFiles):
            myTag = 1
        seq = copy.copy(allSeqFiles[workIDX])
        essc.CreateDirectoryStructure(seq)
        try:
            essc.Run_Modeller(seq,myTag)
        except Exception as err:
            print thisHost,'- Driver error: {0}'.format(err)
            print err, messages[tag],' for index:',workIDX
    def prepFrodan():
        myTag = 0
        seqID = int(workIDX/numModels)
        if seqID >= len(myseqFiles):
            myTag = 1
        seq = copy.copy(allSeqFiles[seqID])
        model = 1 + workIDX%numModels
        try:
            essc.PrepFrodaN(model,seq,myTag)
        except Exception as err:
            print thisHost,'- Driver error: {0}'.format(err)
            print err, messages[tag],' for index:',workIDX
    def frodan():
        myTag = 0
        seqI = (workIDX/(numModels))%len(allSeqFiles)
        seq = allSeqFiles[seqI]
        state = states[workIDX/(len(allSeqFiles) * numModels)]
        if seqI >= len(myseqFiles):
            myTag = 1
        try:
            essc.RunFrodaN(1 + workIDX%numModels,seq,state,myTag)
        except Exception as err:
            print thisHost,'- Driver error: {0}'.format(err)
            print err, messages[tag],' for index:',workIDX
    def prepComplex():
        model = 1 + workIDX%numModels
        Complex = mycomplexStructures[(workIDX/numModels)%len(mycomplexStructures)]
        try:
            essc.PrepComplex(model,Complex)
        except Exception as err:
            print thisHost,'- Driver error: {0}'.format(err)
            print err,':', messages[tag],' for index:',workIDX
    def dockComplex():
        model = 1 + workIDX%numModels
        Complex = mycomplexStructures[(workIDX/numModels)%len(mycomplexStructures)]
        try:
            essc.DockComplex(model,Complex)
        except Exception as err:
            print thisHost,'- Driver error: {0}'.format(err)
            print err,':', messages[tag],' for index:',workIDX
    def mergePDB():
        model = 1 + workIDX%numModels
        Complex = mycomplexStructures[(workIDX/numModels)%len(mycomplexStructures)]
        try:
            essc.MergeAllPDB(str(model).zfill(2),Complex)
        except Exception as err:
            print thisHost,'- Driver error: {0}'.format(err)
            print err,':', messages[tag],' for index:',workIDX
    def genPsize():
        numSequences = len(allStructures)
        seqI = (workIDX/(numModels * pHgridLen*2))%numSequences
        index =  (workIDX%(pHgridLen*2))%pHgridLen
        gridLine = essc.cgridsubsetList[index]
        model = (workIDX/(pHgridLen*2))%numModels
        state = states[(workIDX%(pHgridLen*2))/pHgridLen]
        seq = allStructures[seqI]
        try:
            essc.GenPsize(model,seq,state,gridLine)
        except Exception as err:
            print thisHost,'- Driver error: {0}'.format(err)
            print err, 'for index:',workIDX
    def complexIn():
        numSequences = len(mycomplexStructures)
        seqI = (workIDX/(numModels * pHgridLen))%numSequences
        index =  workIDX%pHgridLen
        Complex = mycomplexStructures[seqI]
        gridLine = essc.cgridsubsetList[index]
        model = (workIDX/pHgridLen)%numModels
        try:
            essc.genComplexIN(model,Complex,gridLine)
        except Exception as err:
            print thisHost,'- Driver error: {0}'.format(err)
            print err, 'for index:',workIDX
    def genIn():
        numSequences = len(allStructures)
        seqI = (workIDX/(numModels * pHgridLen*2))%numSequences
        index =  (workIDX%(pHgridLen*2))%pHgridLen
        gridLine = essc.cgridsubsetList[index]
        model = (workIDX/(pHgridLen*2))%numModels
        state = states[(workIDX%(pHgridLen*2))/pHgridLen]
        seq = allStructures[seqI]
        try:
            essc.genIN(model,seq,state,gridLine)
        except Exception as err:
            print thisHost,'- Driver error: {0}'.format(err)
            print err, 'for index:',workIDX
    def genSAS():
        numSequences = len(allStructures)
        seqI = (workIDX/(numModels * pHgridLen*2))%numSequences
        index =  (workIDX%(pHgridLen*2))%pHgridLen
        gridLine = essc.cgridsubsetList[index]
        model = (workIDX/(pHgridLen * 2))%numModels
        state = states[(workIDX%(pHgridLen*2))/pHgridLen]
        seq = allStructures[seqI]
        try:
            essc.GenSAS(model,seq,state,gridLine)
        except Exception as err:
            print thisHost,'- Driver error: {0}'.format(err)
            print err, 'for index:',workIDX
    def calcPots():
        numSequences = len(allStructures)
        seqI = (workIDX/(numMods * pHgridLen*2))%numSequences
        index =  (workIDX%(pHgridLen*2))%pHgridLen
        gridLine = essc.cgridsubsetList[index]
        model = (workIDX/(pHgridLen*2))%numMods
        prefix,conc,pH = gridLine.split()
        state = states[(workIDX%(pHgridLen*2))/pHgridLen]
        seq = allStructures[seqI]
        try:
            essc.CalcPots(model,seq,state,gridLine)
        except Exception as err:
            print thisHost,'- Driver error: {0}'.format(err)
            print err, 'for index:',workIDX
    def calcBE():
        seq = allSeqFiles[workIDX]
        try:
            essc.Calc_Means_Ratios(seq,numMods)
        except Exception as err:
            print thisHost,'- Driver error: {0}'.format(err)
            print err, 'for index:',workIDX
    def calcMnR():
        Complex = mycomplexStructures[workIDX]
        try:
            essc.Calc_Binding_Energies(Complex,numMods)
        except Exception as err:
            print thisHost,'- Driver error: {0}'.format(err)
            print err, 'for index:',workIDX
    def calcResPots():
        seqI = workIDX/(numMods * pHgridLen *  2)
        index = (workIDX%(pHgridLen*2))%pHgridLen
        gridLine = essc.cgridsubsetList[index]
        model = (workIDX/(pHgridLen*2))%numMods
        state = states[(workIDX%(pHgridLen*2))/pHgridLen]
        seq = allSeqFiles[seqI]
        try:
            essc.CalcRESPots(model,seq,state,gridLine)
        except Exception as err:
            print thisHost,'- Driver error: {0}'.format(err)
            print err, 'for index:',workIDX
    def calResMnR():
        seqI = workIDX/2
        seq = allSeqFiles[seqI]
        state = states[workIDX%2]
        try:
            essc.CalcResMnR(seq,numMods,state)
        except Exception as err:
            print thisHost,'- Driver error: {0}'.format(err)
            print err, 'for index:',workIDX
    def bad():
        print "Bad value handed to 'runStep()'\n"
    switch = {
        0: Complex,
        1: modeller,
        2: prepFrodan,
        3: frodan,
        4: prepComplex,
        5: dockComplex,
        6: mergePDB,
        7: genPsize,
        8: complexIn,
        9: genIn,
        10: genSAS,
        11: calcPots,
        12: calcBE,
        13: calcMnR,
        14: calcResPots,
        15: calResMnR
    }
    def runStep(tag):
        step = switch.get(tag, bad)
        step()
    #### !!!!!!!!!!!!!! ####
    # Begin processing assigned work units
    while tag != 200:
        ESSCrequest = comm.Irecv(MPIbuffer,source=0,tag=MPI.ANY_TAG)
        comm.Send(MPIbuffer,dest=0,tag=10)
        ESSCrequest.Wait(ESSCstatus)
        workIDX = MPIbuffer[0]
        tag = ESSCstatus.Get_tag()
        if tag < 100:
            print '%s for %s' % (messages[tag],workIDX)
            print time.ctime()
            runStep(tag)
        elif tag == 100:
            print 'Finished available work'
            print time.ctime()
            comm.Barrier()
comm.Barrier()

if rank == 0:
    if(int(myconfig['Global']['doRAnalysis'])):
        print 'Performing BESI/EVM Analysis'
        print time.ctime()
        try:
            essc.R_Analysis()
        except Exception as err:
            print thisHost,'- Driver error: {0}'.format(err)
            print err, 'BESI/EVM Analysis'
        print 'Processing complete'
else:
    print 'All assigned work completed'

print time.ctime()

MPI.Finalize()
