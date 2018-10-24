#!/bin/bash

# untar your R installation
tar -xzf R-2018-09-30a.tar.gz
tar -xzf SLIBS.tar.gz
# make sure the script will use your R installation
export PATH=$(pwd)/R/bin:$PATH
export LD_LIBRARY_PATH=$(pwd)/SS:$LD_LIBRARY_PATH
# run R, with the name of your  R script
R CMD BATCH '--args argname='$1' run_num='$2'' Rscript/pvl-chr19-10mb.R 'pvl_run'$2'-'$1'.Rout'


#rm *.Rout
