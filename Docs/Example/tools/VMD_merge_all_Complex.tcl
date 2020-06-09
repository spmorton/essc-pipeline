#==============================================================================
#     VMD_merge_all_Complex.tcl
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
# 10/2/2017
#
# This VMD script merges the supplied list of pdb's (4) into a single pdb 
# of complex structures into one
# 
# command format
# vmd options -e script_file -args args
#
# example command:
# REQUIRED - 1st arg is the absolute path to the Structures folder
# vmd -dispdev none -e VMD_merge_all_Complex.tcl -args mol1 mol2 mol3 mol4

if {$argc < 4} {
    puts "REQUIRED - mol1 mol2 mol3 mol4 the four components of a complex BE calculation using absolute paths"
    q
    }

package require topotools 

set m1 [lindex $argv 0]; list
set m2 [lindex $argv 1]; list
set m3 [lindex $argv 2]; list
set m4 [lindex $argv 3]; list
#set numModels [lindex $argv 1]; list
#cd $wDir
#set dirs [glob -type d Complex-*]; list
#set numDirs [llength $dirs]; list

#set thisPDB "frodan/01/bound/merged.pdb"; list
set mol1 [mol new $m1]; list
set mol2 [mol new $m2]; list
set mol3 [mol new $m3]; list
set mol4 [mol new $m4]; list

#set x [string is integer -strict numModels]; list
#incr x
#format "%02d" 1


#set model [format "%02d" $i]; list
#set nextPDB [concat frodan/$model/bound/merged.pdb]; list
#set mol2 [mol new $nextPDB]; list

# select all atoms from both
set sel1 [atomselect $mol1 all]; list
set sel2 [atomselect $mol2 all]; list
set sel3 [atomselect $mol3 all]; list
set sel4 [atomselect $mol4 all]; list

# Merge the two into one
set tmp1 [::TopoTools::selections2mol "$sel1 $sel2"]; list
set sel5 [atomselect $tmp1 all]; list
set tmp2 [::TopoTools::selections2mol "$sel5 $sel3"]; list
set sel6 [atomselect $tmp2 all]; list
set molOut [::TopoTools::selections2mol "$sel6 $sel4"]; list
#mol delete $mol1
#mol delete $mol2
#set mol1 $mol3; list

animate write pdb "merged.pdb" $molOut

quit
