#! /bin/bash
# Last edited on 2012-02-07 02:50:06 by stolfilocal

funfile="$1"; shift;
dtag="$1"; shift;
stag="$1"; shift;
ymin="$1"; shift;
ymax="$1"; shift;

tmp="/tmp/$$"
outimage="${funfile%%.*}-${dtag}${stag}.png"
tmpimage="${tmp}.png"
  
export GDFONTPATH="."

gnuplot <<EOF
set term png truecolor rounded linewidth 2.0 size 3200,800 font "arial,18"
set output "${tmpimage}"

call "${funfile}"

set yrange [${ymin}:${ymax}]
unset key
set samples 1000

# Max pulse half-width:
HW = 6.0

# Goal function half-width
HG = HW + 2.0

# Clip goal functions to this half-width:
HC = HG + HW + 0.5;

# Spacing between goal function samples:
D = 2*HC + 2.0;

# The unit dirac pulse:
PD(x) = (abs(x) > HC ? 0/0 : (abs(x) > 0.01 ? 0 : 1))

# Polynomials to interpolate and their derivatives, clipped to [-HG:+HG]:
D0P0(x) = (abs(x) > HC ? 0/0 : (abs(x) > HG ? 0 : 1))
D0P1(x) = (abs(x) > HC ? 0/0 : (abs(x) > HG ? 0 : (x/HG)))
D0P2(x) = (abs(x) > HC ? 0/0 : (abs(x) > HG ? 0 : (x/HG)**2))
D0P3(x) = (abs(x) > HC ? 0/0 : (abs(x) > HG ? 0 : (x/HG)**3))

PP(x,r) = (r == 0 ? D0P0(x) : (r == 1 ? D0P1(x) : (r ==2 ? D0P2(x) : (r == 3 ? D0P3(x) : 0/0))))

D1P0(x) = (abs(x) > HC ? 0/0 : 0)
D1P1(x) = (abs(x) > HC ? 0/0 : (abs(x) > HG ? 0 : 1.0/HG))
D1P2(x) = (abs(x) > HC ? 0/0 : (abs(x) > HG ? 0 : 2.0*(x/HG)/HG))
D1P3(x) = (abs(x) > HC ? 0/0 : (abs(x) > HG ? 0 : 3.0*(x/HG)**2/HG))

D2P0(x) = (abs(x) > HC ? 0/0 : 0)
D2P1(x) = (abs(x) > HC ? 0/0 : 0)
D2P2(x) = (abs(x) > HC ? 0/0 : (abs(x) > HG ? 0 : 2.0/HG**2))
D2P3(x) = (abs(x) > HC ? 0/0 : (abs(x) > HG ? 0 : 6.0*(x/HG)/HG**2))

D3P0(x) = (abs(x) > HC ? 0/0 : 0)
D3P1(x) = (abs(x) > HC ? 0/0 : 0)
D3P2(x) = (abs(x) > HC ? 0/0 : 0)
D3P3(x) = (abs(x) > HC ? 0/0 : (abs(x) > HG ? 0 : 6/HG**3))

# Goal functions to interpolate:
G0(x) = ${dtag}P0(x)
G1(x) = ${dtag}P1(x)
G2(x) = ${dtag}P2(x)
G3(x) = ${dtag}P3(x)

GG(x,r) = (r == 0 ? G0(x) : (r == 1 ? G1(x) : (r ==2 ? G2(x) : (r == 3 ? G3(x) : 0/0))))

# The function to plot is {F(x)}:
F(x) = ${dtag}${stag}(x)

# Interpolated polynomials with {F} around 0:
TT(x,r,m) = (m == 0 ? (r == 0 ? PP(0,r)*F(x) : 0) : PP(-m,r)*F(x+m) + PP(+m,r)*F(x-m) + TT(x,r,m-1))

T0(x) = TT(x,0,HG)
T1(x) = TT(x,1,HG)
T2(x) = TT(x,2,HG)
T3(x) = TT(x,3,HG)

set xrange [-HC-2:+4*D+HC+2]

plot \
  (F(x))                                 with lines  lt 1      lc rgb '#000088', \
  "dots.dat" using 1:(PD(column(1)))     with points pt 7 ps 2 lc rgb '#000000', \
  (T0(x-1*D))                            with lines  lt 1      lc rgb '#ff0000', \
  (G0(x-1*D))                            with lines  lt 1      lc rgb '#777777', \
  "dots.dat" using 1:(G0(column(1)-1*D)) with points pt 7 ps 2 lc rgb '#000000', \
  (T1(x-2*D))                            with lines  lt 1      lc rgb '#008800', \
  (G1(x-2*D))                            with lines  lt 1      lc rgb '#777777', \
  "dots.dat" using 1:(G1(column(1)-2*D)) with points pt 7 ps 2 lc rgb '#000000', \
  (T2(x-3*D))                            with lines  lt 1      lc rgb '#ff0000', \
  (G2(x-3*D))                            with lines  lt 1      lc rgb '#777777', \
  "dots.dat" using 1:(G2(column(1)-3*D)) with points pt 7 ps 2 lc rgb '#000000', \
  (T3(x-4*D))                            with lines  lt 1      lc rgb '#008800', \
  (G3(x-4*D))                            with lines  lt 1      lc rgb '#777777', \
  "dots.dat" using 1:(G3(column(1)-4*D)) with points pt 7 ps 2 lc rgb '#000000'
quit
EOF

convert ${tmpimage} -resize '50%' ${outimage}
display ${outimage}
rm -f ${tmpimage}

