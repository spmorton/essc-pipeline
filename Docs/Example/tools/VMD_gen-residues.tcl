#==============================================================================
#     VMD_gen_residues.tcl
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
# 11/23/2017
#
# This VMD script  calculates the SAS and SASA for each residue of a protein
# 
# command format
# vmd options -e script_file -args args
#
# example command:
# vmd -dispdev none -e VMD_gen-residues.tcl -args thisPQRfile.pqr OMPThreads

if {$argc < 3} {
    puts "Example command"
    puts "vmd -dispdev none -e VMD_gen-residues.tcl -args thisPQRfile.pqr OMPThreads"
    return 0
    }

mkdir residues
# set the INPUT pqr file name & number of threads
set filename [lindex $argv 0]
set seq [lindex $argv 1]
set OMPThreads [lindex $argv 2]

# split for the protein/solution designation
set names [split $filename .]
set proteinName [lindex $names 0]

puts "Working on molecule: $seq $filename"

# set the env
set env(VMDFORCECPUCOUNT) $OMPThreads

# load and select the molecule
set mol1 [mol new $filename]
set sel1 [atomselect $mol1 all]

# get the resID list
set resIDS [lsort -integer -unique [$sel1 get resid]]

# loop through the resID list
foreach x $resIDS {
    set res [atomselect $mol1 "resid $x"]
    set myNames [$res get resname]
    set myName [lindex $myNames 0]
    set mysas [measure sasa 1.4 $sel1 -points mygrid -restrict $res]
    #if { [expr double($mysas)] > 0.0 } {
    set gridout ""
    foreach y $mygrid {
        if { [llength $gridout] > 0 } { 
            set gridout "$gridout\n$y"
            } else {
                set gridout "$y"
                }
            }
    set myvmdf "residues/$seq-$proteinName-res-$x-$myName.sas"
    set mysasf "residues/$seq-$proteinName-res-$x-$myName.sasa"
    puts "Writing files: $myvmdf \$mysasf"
    set myf [open $myvmdf w]
    puts $myf $gridout
    close $myf
    set myf [open $mysasf w]
    puts $myf $mysas
    close $myf
        #}
}

quit

