#! /bin/sed -f
# Last edited on 2015-04-29 22:02:44 by stolfilocal

s:# DT_INI_SG:#         COEF  DT_INI_SG:g
s:# ----------:# ------------  ----------:g
s:^ *\(20[01][0-9][-]\):        0.0000  \1:g
