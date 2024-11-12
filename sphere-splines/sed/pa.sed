#! /bin/sed -f
# Last edited on 2005-10-02 22:24:20 by stolfi

# s/SPParams/SPOptions/g
s/SPOptions_T\([ ;{}(),]\)/SPOptions_Parser_t\1/g
s/SPOptions_NewT/SPOptions_NewParser/g
