#==============================================================================
#     VMD_align.tcl
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
# 11/1/2017
#
# This VMD script aligns a protein with another similar protein
# by selecting the first 400 residues of the backbone
# 
# command format
# vmd options -e script_file -args args
#
# example command:
# vmd -dispdev none -e VMD_seperate.tcl -args to_be_fitted.pdb target.pdb output.pdb

if {$argc < 1} {
    puts "REQUIRED - vmd options -e script_file -args args"
    puts "Example command"
    puts "vmd -dispdev none -e VMD_align.tcl -args to_be_fitted.pdb target.pdb output.pdb"
    exit
    }
 
# set the pdb file names
set ml0 [lindex $argv 0]
set ml1 [lindex $argv 1]

set outfile [lindex $argv 2]

set m0 [mol new $ml0]
set m1 [mol new $ml1]

set m0all [atomselect $m0 all]
set sel0 [atomselect $m0 "backbone and resid 0 to 400"]
set sel1 [atomselect $m1 "backbone and resid 0 to 400"]

set M [measure fit $sel0 $sel1]
$m0all move $M

$m0all writepdb $outfile
q
