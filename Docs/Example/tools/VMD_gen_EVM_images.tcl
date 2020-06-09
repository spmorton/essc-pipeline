#==============================================================================
#     VMD_gen_EVM_images.tcl
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
# 4/4/2018
#
# Requires ImageMagick (coded against 6.8.8-1, probably not a hard dependency)
#
# This VMD script  generates all the EVM imagery
# Adjust the rotation matrix and scale settings to your liking by 
# adjusting the position of the protein to your liking and execute
#   molinfo mol# get rotate_matrix
# and pasting the result below
# 
# For targetBase and srcDir these can be the same or different depending on the need
# such that you may have multiple source directories to pull proteins from but want 
# a single destination folder for all output as in the R script that generates the 
# 'Sequence_EVM_Selections.txt' file.

# !!!! ### USER EDITABLE SECTION #### !!!!
# Position a source protein as you want it to be viewed and execute
# molinfo mol# get rotate_matrix
# Paste rotation matrix data here before execution
set rmat "{{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}"
# !!!! ### END USER EDITABLE SECTION #### !!!!

if {$argc < 2} {
    puts "Example command"
    puts "vmd -dispdev none -e VMD_gen_EVM_images.tcl -args explicit_targetBase explicit_srcDir"
    return 0
    }
    
# set the directory base to write images
set targetBase [lindex $argv 0]; list

# set the source directory base
set srcDir [lindex $argv 1]; list

color Display Background white
axes location off

set infile [open "$targetBase/Analysis/Datasets/Sequence_EVM_Selections.txt" r]; list
set Data [read $infile]; list
close $infile 
set DS [split $Data "\n"]; list

set items [llength $DS]; list

for {set i 0} {$i < $items} {incr i} {
    set thisSequence [lindex $DS $i]
    set thisPDB "$srcDir/Structures/$thisSequence/frodan/01/unbound/unbound.pdb"
    set thisMol [mol new $thisPDB]
    set baseFileName "$targetBase/Analysis/EVM_images/$thisSequence"
    set fileName "$baseFileName.tga"
    set fileName2 "$baseFileName.pdf"
    molinfo $thisMol set rotate_matrix $rmat
    scale to 0.025
    mol modstyle 0 $thisMol NewCartoon
    mol modcolor 0 $thisMol structure
    mol addrep $thisMol
    incr i
    mol modselect 1 $thisMol resid [lindex $DS $i]
    mol modstyle 1 $thisMol QuickSurf
    mol modmaterial 1 $thisMol Transparent
    mol modcolor  1 $thisMol ColorID 1
    incr i
    render TachyonInternal $fileName
    # Give the OS a second to write the file out to disk
    sleep 1
    mol off $thisMol
    mol delete $thisMol
    convert $fileName $fileName2
}

quit

