#! /bin/bash
# Last edited on 2014-07-22 19:17:47 by stolfilocal

../gch_point_generator points.txt 50 2 100.0
ls -l points.txt
#cp -av DATA/subspace/quadrado1.txt points.txt
#cp -av DATA/subspace/quadrado2.txt points.txt

../gch_giftwrap out/points.txt out/graph.dot out/faces.txt
# o programa neato pode ser obtido em
#     http://www.research.att.com/sw/tools/graphviz/download.html
neato -Goverlap=false -Tgif out/graph.dot -o out/graph.gif
