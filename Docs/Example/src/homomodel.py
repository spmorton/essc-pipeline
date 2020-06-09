#!/usr/bin/env python

#==============================================================================
#     homomodel.py - customized version for ESSC pipeline
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


# Written for the ESSC Pipeline library
# from a copy provided by modeller
# by Scott P Morton 2020

import sys, os, pickle
from glob import glob
from modeller import *
from modeller.automodel import *

log.level(output=0, notes=0, warnings=0, errors=1, memory=0)
log.none()
#log.verbose()
#log.minimal()
env = environ()
env.io.atom_files_directory = ['./','../atom_files/']
aln = alignment(env)

## If you have multiple chains, multiple templates...
files=glob('template?.pdb')
codes=[]
for myfile in (files):
    print (myfile)
    code = myfile.rstrip('.pdb')
    #print('\\n\\n\\n',code,'\\n\\n\\n')
    codes.append(code)
    mdl = model(env, file=myfile)
    aln.append_model(mdl, atom_files=code, align_codes=code)
    
for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True),
                                    ((1., 0.5, 1., 1., 1., 0.), False, True),
                                    ((1., 1., 1., 1., 1., 0.), True, False)):
    aln.salign(rms_cutoff=3.5, normalize_pp_scores=False,
                rr_file='$(LIB)/as1.sim.mat', overhang=30,
                gap_penalties_1d=(-450, -50),
                gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0,
                dendrogram_file='fm00495.tree',
                alignment_type='tree', # If 'progresive', the tree is not
                                        # computed and all structues will be
                                        # aligned sequentially to the first
                feature_weights=weights, # For a multiple sequence alignment only
                                        # the first feature needs to be non-zero
                improve_alignment=True, fit=True, write_fit=write_fit,
                write_whole_pdb=whole, output='ALIGNMENT QUALITY',
                edit_file_ext=('.pdb','_fit.pdb'))

aln.write(file='template.pap', alignment_format='PAP')
aln.write(file='template.ali', alignment_format='PIR')

aln.salign(rms_cutoff=1.0, normalize_pp_scores=False,
            rr_file='$(LIB)/as1.sim.mat', overhang=30,
            gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3),
            gap_gap_score=0, gap_residue_score=0, dendrogram_file='1is3A.tree',
            alignment_type='progressive', feature_weights=[0]*6,
            improve_alignment=False, fit=False, write_fit=True,
            write_whole_pdb=False, output='QUALITY')

## Restart now

## Restart now
env = environ()
env.libs.topology.read(file='$(LIB)/top_heav.lib')

# Read aligned structure(s):
aln = alignment(env)
aln.append(file='template.ali', align_codes='ALL')
aln_block = len(aln)

# Read aligned sequence(s):
aln.append(file='sequence.ali', align_codes='ALL')

# Structure sensitive variable gap penalty sequence-sequence alignment:
aln.salign(output='', max_gap_length=20,
            gap_function=True,   # to use structure-dependent gap penalty
            alignment_type='PAIRWISE', align_block=aln_block,
            feature_weights=(1., 0., 0., 0., 0., 0.), overhang=0,
            gap_penalties_1d=(-450, 0),
            gap_penalties_2d=(0.35, 1.2, 0.9, 1.2, 0.6, 8.6, 1.2, 0., 0.),
            similarity_flag=True)

aln.write(file='model.ali', alignment_format='PIR')
aln.write(file='model.pap', alignment_format='PAP')

a = automodel(env, alnfile = 'model.ali', knowns = codes,
              sequence = 'model', assess_methods=(assess.DOPE, assess.GA341))
a.starting_model= 1
a.ending_model  = 3 * ++!!NUM_MODELS!!++
a.make()

# Get a list of all successfully built models from a.outputs
ok_models = [x for x in a.outputs if x['failure'] is None]

of = open('modeller.log','w')
for x in a.outputs:
    of.write('Model Name - '+str(x['name'])+'\n')
    of.write('Dope score - '+str(x['DOPE score'])+'\n')
    of.write('GA341 score - '+str(x['GA341 score'])+'\n')
    of.write('failure - '+str(x['failure'])+'\n')
    of.write('Model num - '+str(x['num'])+'\n')
    of.write('molpdf - '+str(x['molpdf'])+'\n')
    of.write('pdfterms - '+str(x['pdfterms'])+'\n\n')
of.close()

# Rank the models by DOPE score
key = 'DOPE score'
if sys.version_info[:2] == (2,3):
    # Python 2.3's sort doesn't have a 'key' argument
    ok_models.sort(lambda a,b: cmp(a[key], b[key]))
else:
    ok_models.sort(key=lambda a: a[key])

with open ('models_by_DOPE.pickle','wb') as f:
    pickle.dump(ok_models,f)

for i in range(1,++!!NUM_MODELS!!++ + 1):
    destination = 'ESSC.%s.pdb' % (str(i).zfill(4))
    try:
        os.symlink(ok_models[i-1]['name'],destination)
    except Exception as err:
        print('(homomodel os.symlink) OS error: {0}'.format(err))
        print 'forcing replacement of symlink')
        os.remove(destination)
        os.symlink(ok_models[i-1]['name'],destination)

# Get top modell
m = ok_models[0]
print("Top model: %s (DOPE score %.3f)" % (m['name'], m[key]))
