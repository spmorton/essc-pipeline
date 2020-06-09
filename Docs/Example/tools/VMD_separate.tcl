#==============================================================================
#     VMD_separate.tcl
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
# 4/23/2017
#
# This VMD script opens a sequence specified by command line argument(CLA) 0
# Then selects the gp120 chain and makes that chain the top of the system
# The mode is selected by argument 2 is either 0 or 1 and determines the 
# operation separate or rejoin respectively
# 
# command format
# vmd options -e script_file -args args
#
# example command:
# vmd -dispdev none -e VMD_separate.tcl -args 2ny7.pdb 0

# deprecated but retained for future reference
# # vmd -dispdev none -e VMD_separate.tcl -args 2ny7.pdb 10 10 10

if {$argc < 1} {
    puts "REQUIRED - vmd options -e script_file -args args"
    puts "Example command"
    puts "vmd -dispdev none -e VMD_separate.tcl -args 2ny7.pdb 0"
    exit
    }
 
# set the pdb file name
set seq [lindex $argv 0]; list
# set the mode of operation IE separate(0) or rejoin(1)
set mode [lindex $argv 1]; list

# deprecated
# set the distance for 'moveby' from the last 3 CLA's
# set v "[lindex $argv 1] [lindex $argv 2] [lindex $argv 3]"

# Open the sequence as a new molecule
mol new $seq

# select the entire set of structures
set selall [atomselect top all]; list

# Determine the chains that are present
set chains [lsort -ascii -unique [$selall get chain]]; list

# Remove L and H from the list
set chains [lsearch -all -inline -not -exact $chains L]; list
set chains [lsearch -all -inline -not -exact $chains H]; list

# select the 'G' chain and make it the top of the system
#set sel [atomselect top "chain G"]

# select the gp120 chain and the antibody H chain
set sel [atomselect top "chain $chains"]; list
set selH [atomselect top "chain H"]; list

# Get the locations of each
set loc [measure center $sel weight mass]; list
set locH [measure center $selH weight mass]; list

# Get the distance between the two centers
# basically we will double this or cut in half
set diff [vecsub $loc $locH]; list

# If in rejoin mode cut the distance by .5 and change the sign
# set the file name appropriately for the operation
if {$mode == 1} {
    set diff [vecscale -.5 $diff]; list
    # because we are working with the separated file
    # we need to trim off the format designation
    set spl [split $seq _]; list
    set newfile [lindex $spl 0]; list
    append newfile "_rejoined.pdb"

    } else {
    # split the seq file name and generate the new file name
    set spl [split $seq .]; list
    set newfile [lindex $spl 0]; list
    append newfile "_separated.pdb"
}

# move the selected chain by the given distances
$sel moveby $diff

# write the new PDB file
$selall writepdb $newfile

quit

