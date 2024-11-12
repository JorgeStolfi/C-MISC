#! /bin/bash 
# Last edited on 2004-12-26 11:42:14 by stolfi

this=$$
( sleep 5 && kill -s 15 ${this} ) &

display ~/misc/egino/apes.jpg
kill -s 15 %sleep

