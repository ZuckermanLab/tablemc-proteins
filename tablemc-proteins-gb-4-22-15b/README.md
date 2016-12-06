Tabulated Monte Carlo for proteins, using CHARMM 19 force field and FACTS implicit solvent. (This code is more advanced than the one used for Spiriti and Zuckerman, JCTC 10, 5161 (2014)).  At the moment the table functionality and the "peratomconv" and "perfrag" GB modes (see documentation) are not to be used.

# src/

Contains source files, header files and compilation script.

# data/

Contains force field file, amino acid definitions file, and fragment reference geometries.

## data/charmm19.prm

A Tinker parameter file for the CHARMM 19 force field.

## data/defs-protein5-charmm.txt

"Definitions" file describing how to construct each amino acid from fragments.

## data/fragments3/*.xyz

Fragment reference geometries (one file for each fragment type).

# tablemc-proteins-6-5-15.tex

Instructions on compiling and using the code.

# tablemc-proteins-internals-6-1-15.tex

Brief description of internal operation of the code.


