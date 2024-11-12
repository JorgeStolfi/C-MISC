/* Last edited on 2024-11-12 17:20:57 by stolfi */
/* Vigenère cypher example for wikipedia. */

#include <stdio.h>
#include <string.h>
int main(int argc, char**argv){ 
  fprintf(stderr, "----------------------------------------------------------------------\n");
  char *arg = "attacking tonight"; fprintf(stderr, "%s\n", arg);
  char *key = "oculorhinolaryngology"; fprintf(stderr, "%s\n", key);
  int n = strlen(arg);
  for (int i = 0; i < n; i++)
    { char enc;
      if ((arg[i] >= 'a') && (arg[i] <= 'z'))
        { enc = (((arg[i]-'a') + (key[i]-'a')) % 26) + 'a'; }
      else
        { enc = arg[i]; }
      fputc(enc, stderr);
    }
  fputc('\n', stderr);
  fprintf(stderr, "----------------------------------------------------------------------\n");
  return 0;
}
