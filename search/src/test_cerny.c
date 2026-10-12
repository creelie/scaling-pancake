#include <stdio.h>
#include "rt.h"
int main(){ for(int n=2;n<=20;n++){ DFA A; A.n=n; A.k=2; for(int q=0;q<n;q++){ A.d[0][q]=(q==0)?1:q; A.d[1][q]=(q+1)%n; } dfa_prep(&A); printf("n=%d rt=%d (n-1)^2=%d\n",n,dfa_rt(&A),(n-1)*(n-1)); } }
