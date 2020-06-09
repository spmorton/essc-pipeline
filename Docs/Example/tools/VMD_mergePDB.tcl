#==============================================================================
#     VMD_mergePDB.tcl
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



# Written by Scott P. Morton
# MTSU College of Graduate Studies
# Computational Sciences Program
# 5/16/2017
#
# This VMD script merges two pdb files into one
# 
# command format
# vmd options -e script_file -args args
#
# example command:
# REQUIRED - 1st arg is the gp120 pdb 2cnd is the antibody or other
# vmd -dispdev none -e VMD_mergePDB.tcl -args 2ny7.pdb bnab.pdb

if {$argc < 3} {
    puts "REQUIRED - 1st arg is the gp120 pdb 2cnd is the antibody or other"
    puts "optional outfile name is the 3rd arg"
    puts "pdb's must be pre-aligned through salign or other method"
    puts "Example command"
    puts "vmd -dispdev none -e VMD_mergePDB.tcl -args 2ny7.pdb bnab.pdb myout.pdb"
    q
    }

package require topotools 

# Get the pdb file names from the CLA
set p1 [lindex $argv 0]; list
set p2 [lindex $argv 1]; list

if {$argc == 4} {
    set outfile [lindex $argv 2]; list
    set pdbfile $outfile; list
    } else {
    set pdbfile "merged.pdb"; list
    }

# Load the pdb's
set mol1 [mol new $p1]; list
set mol2 [mol new $p2]; list

# select all atoms from both
set sel1 [atomselect $mol1 all]; list
set sel2 [atomselect $mol2 all]; list

# Merge the two into one
set mol3 [::TopoTools::selections2mol "$sel1 $sel2"]; list

# Write the data structures to disk
# animate write psf $psffile $mol3 
animate write pdb $pdbfile $mol3

quit
