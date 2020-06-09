#==============================================================================
#     VMD_fixchainids.tcl
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
# This VMD script  fixes the chain ID mess that editconf leaves
# 
# command format
# vmd options -e script_file -args args
#
# example command:
# vmd -dispdev none -e VMD_fixchainids.tcl -args merged.pdb merged2.pdb 3 G H L

if {$argc < 5} {
    puts "Example command"
    puts "vmd -dispdev none -e VMD_fixchainids.tcl -args merged.pdb merged2.pdb 3 G H L"
    return 0
    }
    
# set the INPUT pdb file name
set mol1 [lindex $argv 0]; list

# set the OUTPUT pdb file name
set pdbOUT [lindex $argv 1]; list

set numSegs [lindex $argv 2]; list

# Open the sequence as a new molecule
mol new $mol1

# select the entire set of structures
set selall [atomselect top all]; list

if {$numSegs == 3} {
    set seg1 [lindex $argv 3]; list
    set seg2 [lindex $argv 4]; list
    set seg3 [lindex $argv 5]; list

    # select the fragments
    set frg0 [atomselect top "fragment 0"]; list
    set frg1 [atomselect top "fragment 1"]; list
    set frg2 [atomselect top "fragment 2"]; list

    # Correct the chain ID's
    $frg0 set chain $seg1
    $frg1 set chain $seg2
    $frg2 set chain $seg3
    
    } elseif {$numSegs == 2} {
    
    set seg1 [lindex $argv 3]; list
    set seg2 [lindex $argv 4]; list

    # select the fragments
    set frg0 [atomselect top "fragment 0"]; list
    set frg1 [atomselect top "fragment 1"]; list

    # Correct the chain ID's
    $frg0 set chain $seg1
    $frg1 set chain $seg2
    } elseif {$numSegs == 1} {
    
    set seg1 [lindex $argv 3]; list

    # select the fragments
    set frg0 [atomselect top "fragment 0"]; list

    # Correct the chain ID's
    $frg0 set chain $seg1
}    
    
# Write the new PDB file
puts "Writing PDB file $pdbOUT"
$selall writepdb $pdbOUT
quit
