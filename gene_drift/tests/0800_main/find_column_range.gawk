#! /usr/bin/gawk -f
# Last edited on 2023-05-31 20:54:58 by stolfi

# Finds the range of years in a per-range file

# Reads from {stdin} a file. Looks at a specified column (default column 1),
# which should contain numbers, possibly negative.  Writes the min and max 

BEGIN { 
  yMin = +1.0e200; 
  yMax = -1.0e200; 
  if (col == "") { col = 1; }
}

(match($(col), /^[-+]?[0-9]*([0-9]|([0-9][.]|[.][0-9])[0-9]*)(|[Ee][-+]?[0-9]+)$/)) { 
  y = $(col) + 0;
  if (y < yMin) { yMin = y; }
  if (y > yMax) { yMax = y; }
}

END { printf "%d %d\n", yMin, yMax; }


