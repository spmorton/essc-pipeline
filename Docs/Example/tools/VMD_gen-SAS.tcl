#==============================================================================
#     VMD_gen_SAS.tcl
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
# 1/6/2018
#
# This VMD script  calculates the SAS and SASA for a protein
# 
# command format
# vmd options -e script_file -args args
#
# example command:
# vmd -dispdev none -e VMD_gen-SAS.tcl -args thisPQRfile.pqr OMPThreads

if {$argc < 3} {
    puts "Example command"
    puts "vmd -dispdev none -e VMD_gen-SAS.tcl -args thisPQRfile.pqr OMPThreads"
    return 0
    }

# set the INPUT pqr file name & number of threads
set filename [lindex $argv 0]; list
set OMPThreads [lindex $argv 1]; list

# set the env
set env(VMDFORCECPUCOUNT) $OMPThreads; list

# load and select the molecule
set mol1 [mol new $filename]; list
set sel1 [atomselect $mol1 all]; list
set mysas [measure sasa 1.4 $sel1 -points mygrid]; list

# format the grid data
foreach y $mygrid {
        append gridout \n $y
}
set mygrid [ string trim $gridout ]; list

set myvmdf "$filename.sas"; list
set mysasf "$filename.sasa"; list

puts "Working in folder: [pwd] \n$filename"
puts "SASA calculated: $mysas"
puts "Writing files: \n$myvmdf \n$mysasf"
set myf [open $myvmdf w]; list
puts $myf $mygrid
close $myf
set myf [open $mysasf w]; list
puts $myf $mysas
close $myf

quit
