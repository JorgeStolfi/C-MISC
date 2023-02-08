#!/bin/csh -f

###############################
#  Degenerate point files     #
###############################
#cp ../points/subspace/cubo.txt points.txt; set scale=0.6; set distortion=1
#cp ../points/subspace/cubo2.txt points.txt; set scale=0.1; set distortion=1
cp ../points/subspace/cubo3.txt points.txt; set scale=0.1; set distortion=1

./giftwrap points.txt graph.dot faces.txt
# o programa neato pode ser obtido em
#     http://www.research.att.com/sw/tools/graphviz/download.html
~/bin/neato -Goverlap=false -Tgif graph.dot -o graph.gif
./point_proj_hc points.txt pointsProj.txt ${distortion}
./pov_generator_hc pointsProj.txt faces.txt main.pov ${scale}
make -f Makefile_povray
