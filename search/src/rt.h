// Reset threshold of a DFA with n<=28 states, k letters, via BFS over subsets (power automaton).
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
typedef struct { int n,k; int d[8][32]; uint32_t tab[8][4][256]; } DFA;
static void dfa_prep(DFA*A){
  for(int a=0;a<A->k;a++) for(int byte=0;byte<4;byte++) for(int v=0;v<256;v++){
    uint32_t img=0; for(int b=0;b<8;b++) if(v>>b&1){int q=byte*8+b; if(q<A->n) img|=1u<<A->d[a][q];}
    A->tab[a][byte][v]=img; }
}
static inline uint32_t dfa_img(const DFA*A,int a,uint32_t S){
  return A->tab[a][0][S&255]|A->tab[a][1][(S>>8)&255]|A->tab[a][2][(S>>16)&255]|A->tab[a][3][S>>24];
}
// returns reset threshold or -1 if not synchronizing. Uses caller-provided visited array of size 2^n bytes (bitset ok)
static uint8_t *VIS=NULL; static uint32_t *QU=NULL; static int VISN=-1;
static int dfa_rt(const DFA*A){
  int n=A->n; size_t N=(size_t)1<<n;
  if(VISN!=n){ free(VIS); free(QU); VIS=calloc(N/8+1,1); QU=malloc(sizeof(uint32_t)*N); VISN=n; }
  else memset(VIS,0,N/8+1);
  uint32_t full=(n==32)?0xffffffffu:((1u<<n)-1);
  size_t h=0,t=0; QU[t++]=full; VIS[full>>3]|=1<<(full&7);
  int depth=0; size_t levelEnd=t;
  if(__builtin_popcount(full)==1) return 0;
  while(h<t){
    uint32_t S=QU[h++];
    for(int a=0;a<A->k;a++){ uint32_t T=dfa_img(A,a,S);
      if(!(VIS[T>>3]>>(T&7)&1)){ VIS[T>>3]|=1<<(T&7); if(__builtin_popcount(T)==1) return depth+1; QU[t++]=T; } }
    if(h==levelEnd){ depth++; levelEnd=t; }
  }
  return -1;
}
