/* Last edited on 2012-12-15 23:23:11 by stolfilocal */
/* Checks the multiplier circuits. */

#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <sign.h>

#include <logsys.h>
#include <logsys_arith.h>

int main(int argc, char **argv)
  {

    int nx = 12;      /* Number of bits in {X}. */
    int ny = 24;      /* Number of bits in {Y}. */
    int nz = nx + ny; /* Number of bits in {Z=X*Y}. */

    /* Since we use 32/64 bit integers: */
    assert(nx <= 32);
    assert(ny <= 32);

    /* Create the basic system: */
    logsys_va_t *xb[nx]; /* {xb[0..nx-1]} are the bits of {X}. */
    logsys_va_t *yb[ny]; /* {yb[0..ny-1]} are the bits of {Y}. */
    logsys_va_t *zb[nz]; /* {zb[0..nz-1]} are the bits of {Z}. */
    logsys_t *S = logsys_build_mul(nx, ny, xb, yb, zb);
    logsys_print_system(stdout, "--- multiplier ---", S, "----------------");
    logsys_check_mul(S, nx, ny, xb, yb, zb);

    return 0;
  }
