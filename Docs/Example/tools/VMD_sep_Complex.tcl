#==============================================================================
#     VMD_sep_Complex.tcl
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
# 5/18/2017
#
# This VMD script seperates the Env from the 2 chain Abs
#   chains must be G H L
# 
# command format
# vmd options -e script_file -args args
#
# example command:
# REQUIRED - source_pdb target_pdb target_pdb
# vmd -dispdev none -e VMD_sep_Complex.tcl -args 2ny7.pdb /test/2ny7/01/bound/bound.pdb /test/b12/01/bound/bound.pdb

if {$argc < 4} {
    puts "REQUIRED -"
    puts "1st arg is the complex pdb to split"
    puts "2cnd arg is type of secondary component AB or CD4"
    puts "3rd is the absolute path to store the Env"
    puts "4th is the absolute path to store the AB"
    puts "Example command"
    puts "vmd -dispdev none -e VMD_sep_Complex.tcl -args 2ny7.pdb CD4 /test/2ny7/01/bound/bound.pdb /test/b12/01/bound/bound.pdb"
    q
    }

package require topotools 

# Get the pdb file names from the CLA
set p1 [lindex $argv 0]; list
set stype [lindex $argv 1]; list
set Envfile [lindex $argv 2]; list
set ABfile [lindex $argv 3]; list


# Load the pdb's
set mol1 [mol new $p1]; list

# select all atoms from both
set sel1 [atomselect $mol1 "chain G"]; list
if {$stype == 'AB'} {
    set sel2 [atomselect $mol1 "chain H L"]; list
    } elseif { $stype == 'CD4'} {
    set sel2 [atomselect $mol1 "chain C"]; list
    }
    
# Write the data structures to disk
# animate write psf $psffile $mol3 
$sel1 writepdb $Envfile
#animate write pdb $ABfile $ABmol
$sel2 writepdb $ABfile

quit
