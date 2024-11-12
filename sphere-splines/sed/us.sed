#! /bin/sed -f
# Last edited on 2005-10-02 20:47:13 by stolfi

/[#]define PPUSAGE/d
s/PPUSAGE/SPParams_SetUsage/g
