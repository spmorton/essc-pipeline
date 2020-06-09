#!/bin/bash

if [ "$#" != "2" ]; then
    echo
    echo "Usage: $0 [reference.pdb] [fitting.pdb]"
    echo "   Uses sequence and structure information to"
    echo "   produce and optimal alignment between two"
    echo "   structures."
    echo
    exit -1
fi

PDB1="$1"
PDB2="$2"

cat<<EOF | modpy.sh python
# Illustrates the SALIGN multiple structure/sequence alignment

from modeller import *

log.verbose()
env = environ()
env.io.atom_files_directory = ['./','../atom_files/']

aln = alignment(env)
codes=[]

## If you have multiple chains, multiple templates...
for (code) in ('${PDB1}'.rstrip('.pdb'),'${PDB2}'.rstrip('.pdb')):
    codes.append(code)
    mdl = model(env, file=code+'.pdb')
    aln.append_model(mdl, atom_files=code, align_codes=code)

for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True),
                                    ((1., 0.5, 1., 1., 1., 0.), False, True),
                                    ((1., 1., 1., 1., 1., 0.), True, False)):
    aln.salign(rms_cutoff=3.5, normalize_pp_scores=False,
               rr_file='\$(LIB)/as1.sim.mat', overhang=30,
               gap_penalties_1d=(-450, -50),
               gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0,
               dendrogram_file=codes[0]+'_'+codes[1]+'_fm00495.tree',
               alignment_type='tree', # If 'progresive', the tree is not
                                      # computed and all structues will be
                                      # aligned sequentially to the first
               feature_weights=weights, # For a multiple sequence alignment only
                                        # the first feature needs to be non-zero
               improve_alignment=True, fit=True, write_fit=write_fit,
               write_whole_pdb=whole, output='ALIGNMENT QUALITY',
               edit_file_ext=('.pdb','_fit.pdb'))

aln.write(file=codes[0]+'_'+codes[1]+'.pap', alignment_format='PAP')
aln.write(file=codes[0]+'_'+codes[1]+'.ali', alignment_format='PIR')

aln.salign(rms_cutoff=1.0, normalize_pp_scores=False,
           rr_file='\$(LIB)/as1.sim.mat', overhang=30,
           gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3),
           gap_gap_score=0, gap_residue_score=0, 
           dendrogram_file=codes[0]+'_'+codes[1]+'_1is3A.tree',
           alignment_type='progressive', feature_weights=[0]*6,
           improve_alignment=False, fit=False, write_fit=True,
           write_whole_pdb=False, output='QUALITY')

EOF
