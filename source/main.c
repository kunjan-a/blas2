#include<stdio.h>
#include"blas.h"

int main()
{ int i,j;
  float A[5][5];
  float **p;
  p=A;
  float X[5],Y[5];
  float ALPHA,BETA;
  float result[5];
  for(i=0;i<5;i++)
   {
   for(j=0;j<5;j++)
      {A[i][j]=j;
      }
   X[i]=i;
   Y[i]=i;
   }
   ALPHA=1;
   BETA=1;
   result =SGEMV('N',5,5,ALPHA,p,0,X,1,BETA,Y,1);
   for(i=0;i<5;i++)
     printf(" %lf ",result[i]);

}        
  
