#!/bin/bash -l
#source /etc/profile.d/modules.sh
eval `modulecmd bash del python`
eval `modulecmd bash add modeller gromacs vmd apbs frodan pdb2pqr python-2.6 raxml mafft`
eval `modulecmd bash del gromacs`
eval `modulecmd bash add gromacs-5.1.2`
ulimit -u 4096
ulimit -n 4096
python $@

