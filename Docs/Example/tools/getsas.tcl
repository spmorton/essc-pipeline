
#==============================================================================
#     getsas.tcl
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

set files prot-3.8.pqr

if {$files != ""} {
    foreach filename $files {
        puts "Working on molecule: $filename"
        set env(VMDFORCECPUCOUNT) '8'
        set mymol [mol new $filename]
        set mysel [atomselect $mymol all]
        set mysas [measure sasa 1.4 $mysel -points mygrid]
        puts "SASA calculated: $mysas"
        set myvmdf "$filename.vmd"
        set mysasf "$filename.sasa"
        puts "Writing files: $myvmdf \$mysasf"
        set myf [open $myvmdf w]
        puts $myf $mygrid
        close $myf
        set myf [open $mysasf w]
        puts $myf $mysas
        close $myf
    }
}
quit

