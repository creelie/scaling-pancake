// Simulated annealing over k-letter n-state DFAs maximizing reset threshold.
// Reports any automaton with rt > (n-1)^2 (a counterexample to Cerny) and the best non-Cerny-isomorphic score.
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "rt.h"
static uint64_t rs=88172645463325252ull; static inline uint64_t rnd(){rs^=rs<<13;rs^=rs>>7;rs^=rs<<17;return rs;}
static double urand(){return (rnd()>>11)*(1.0/9007199254740992.0);}
int main(int argc,char**argv){
  int n=atoi(argv[1]),k=atoi(argv[2]); long iters=atol(argv[3]); rs^=(uint64_t)atol(argv[4])*0x9E3779B97F4A7C15ull; double T0=argc>5?atof(argv[5]):3.0;
  int bound=(n-1)*(n-1); DFA A; A.n=n;A.k=k; int best=-1; long restarts=0;
  for(;;){
    // random start
    for(int a=0;a<k;a++)for(int q=0;q<n;q++)A.d[a][q]=rnd()%n; dfa_prep(&A); int cur=dfa_rt(&A); if(cur<0)cur=0;
    for(long it=0;it<iters;it++){
      double T=T0*(1.0-(double)it/iters)+0.05;
      int a=rnd()%k,q=rnd()%n,old=A.d[a][q],nw=rnd()%n; if(nw==old)continue;
      A.d[a][q]=nw; dfa_prep(&A); int r=dfa_rt(&A); int sc=r<0?0:r;
      if(sc>=cur || urand()<exp((sc-cur)/T)){ cur=sc;
        if(sc>best){best=sc; printf("n=%d k=%d rt=%d bound=%d restart=%ld it=%ld :",n,k,sc,bound,restarts,it);
          for(int b=0;b<k;b++){printf(" [");for(int p=0;p<n;p++)printf("%d%s",A.d[b][p],p<n-1?",":"");printf("]");} printf("\n"); fflush(stdout);
          if(sc>bound){printf("COUNTEREXAMPLE\n");fflush(stdout);} }
      } else {A.d[a][q]=old;}
    }
    restarts++;
  }
}
