/* See {csc_stock_sheet.h}. */
/* Last edited on 2020-01-01 02:58:07 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <values.h>
#include <string.h>
#include <assert.h>

#include <vec.h>
#include <argparser.h>
#include <r2.h>
#include <bool.h>
#include <affirm.h>

#include <csc_stock_sheet.h>

vec_typeimpl(csc_stock_sheet_vec_t, csc_stock_sheet_vec, csc_stock_sheet_t);

csc_stock_sheet_vec_t csc_stock_sheet_parse_options(argparser_t *pp)
  {
    int32_t ns = 0; /* Number of stock sheets, including repetitions. */
    csc_stock_sheet_vec_t shv = csc_stock_sheet_vec_new(10);
    while (argparser_keyword_present(pp, "-stockSheets"))
      { char *tag = argparser_get_next_non_keyword(pp);
        int32_t nr = (int32_t)argparser_get_next_int(pp, 1, 100);
        char *mat = argparser_get_next_non_keyword(pp);
        double thk = argparser_get_next_double(pp, 0.001, 100.0);
        r2_t size;
        size.c[0] = argparser_get_next_double(pp, 10.0, 10000.0);
        size.c[1] = argparser_get_next_double(pp, 10.0, 10000.0);
        for (uint32_t k = 0;  k < nr; k++)
          { csc_stock_sheet_vec_expand(&(shv), ns); 
            char *tagk = NULL;
            char *tagk = jsprintf("%s.%d", tag, k+1);
            shv.e[ns] = (csc_stock_sheet_t){ .tag = tagk, .size = size, .mat = mat, .thk = thk };
            ns++;
          }
      }
    csc_stock_sheet_vec_trim(&(shv), ns);
    return shv;
  }
