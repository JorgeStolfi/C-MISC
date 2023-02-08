#!/bin/bash
# Last edited on 2015-12-01 14:29:28 by stolfilocal

###############################
#  3-dimensional point files  #
###############################

#../gch_point_generator_sphere out/points.txt 500 3 1000; scale=0.01
#../gch_point_generator out/points.txt 50 3 100; scale=0.1
#cp -av DATA/d3/simplexo.txt out/points.txt; scale=0.6
#cp -av DATA/d3/cubo.txt out/points.txt; scale=0.6
#cp -av DATA/d3/cubo2.txt out/points.txt; scale=0.1
#cp -av DATA/d3/octaedro.txt out/points.txt; scale=1.0
#cp -av DATA/d3/octaedro2.txt out/points.txt; scale=0.2
cp -av DATA/d3/12-cell.txt out/points.txt; scale=1.0
#cp -av DATA/d3/12-cell2.txt out/points.txt; scale=0.2

#cp -av DATA/tmp/tmp.txt out/points.txt; scale=0.6

#cp -av DATA/d3/cubo_lifted.txt out/points.txt; scale=0.6

../gch_giftwrap out/points.txt out/graph.dot out/faces.txt
# o programa neato pode ser obtido em
#     http://www.research.att.com/sw/tools/graphviz/download.html
neato -Goverlap=false -Tgif out/graph.dot -o out/graph.gif
../gch_pov_generator out/points.txt out/faces.txt povray/main.pov ${scale} 0
( cd povray && make )

