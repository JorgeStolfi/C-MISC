#!/bin/bash
# Last edited on 2014-07-22 12:59:55 by stolfilocal

###############################
#  4-dimensional point files  #
###############################

#../gch_point_generator_sphere out/points.txt 100 4 100; scale=0.1; distortion=1
#../gch_point_generator out/points.txt 40 4 15.0; scale=1.0; distortion=5
#cp -av DATA/d4/simplex.txt out/points.txt; scale=0.7; distortion=5
#cp -av DATA/d4/hipercubo.txt out/points.txt; scale=0.5; distortion=40
#cp -av DATA/d4/hipercubo2.txt out/points.txt; scale=0.10; distortion=6
#cp -av DATA/d4/hiperoctaedro.txt out/points.txt; scale=1.0; distortion=5
cp -av DATA/d4/hiperoctaedro2.txt out/points.txt; scale=0.3; distortion=5
#cp -av DATA/d4/24-cell.txt out/points.txt; scale=1.2; distortion=30
#cp -av DATA/d4/24-cell2.txt out/points.txt; scale=0.25; distortion=5

../gch_giftwrap out/points.txt out/graph.dot out/faces.txt
# o programa neato pode ser obtido em
#     http://www.research.att.com/sw/tools/graphviz/download.html
neato -Goverlap=false -Tgif out/graph.dot -o out/graph.gif
../gch_point_proj_hc out/points.txt out/pointsProj.txt ${distortion}
../gch_pov_generator out/pointsProj.txt out/faces.txt povray/main.pov ${scale} 1
( cd povray && make )
