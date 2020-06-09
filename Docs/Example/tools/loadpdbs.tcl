 
set files [glob ~/PNAS/ptree/*.pdb]

foreach item $files {
    mol new $item
}
