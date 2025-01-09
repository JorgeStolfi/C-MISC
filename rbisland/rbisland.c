/* Last edited on 2024-12-21 11:27:11 by stolfi  */

#include <stdio.h>
#include <stdlib.h>

#include <jsrandom.h>

int main(int k, char**a) {
  int N=atoi(a[1]),G=atoi(a[2]),F=atoi(a[3]),T=atoi(a[4]);
  int c[N],p[N+1],t=T,g,s,i;
  i=N+1;for(;i;)p[--i]=0;
  while(t){
    g=G;
    while(g){
      s=0;
      /* printf("%6d ",g); */
      i=N;while(i--)s+=(c[i]=(g==G?(rand()%N<=F):c[(rand()>>2)%N]));
      /* i=0;for(;i<N;)putchar(48+c[i++]);putchar(10); */
      printf("0 %d %d\n",G-g,s);
      g--;
    }
    /* printf("%6d %6d\n",t,s); */
    putchar(10);
    p[s]++;
    t--;
  }
  i=0;while(i<=N)printf("1 %5d %.6f\n",i,(p[i++]+.0)/T);
  return 0;
}
