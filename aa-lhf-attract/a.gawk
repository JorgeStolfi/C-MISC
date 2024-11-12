BEGIN {
 a=1
 o=1
}

{
 A=2^$1
 a/=2
 print $0,$3/$2, $2/A,$3/A, $3*a,$3/o
 o=$3
}
