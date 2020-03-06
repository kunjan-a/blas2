#include"blas.h"

void SGEMV(char TRANS,int M,int N,float ALPHA,float* A,int LDA,float* X,int INCX,float BETA,float* Y,int INCY)
{ 
	int LENX,LENY;
	int KX,KY,IY,IX,JX,JY;
	int i,j;
	float TEMP;
	if(TRANS=='N'||TRANS=='n') 
        {  LENX = N;
           LENY = M;
        }
    else
        {  LENX = M;
           LENY = N;
        }
	if (INCX>0) 
          KX = 1;
    else
          KX = 1 - (LENX-1)*INCX;
    if (INCY>0) 
          KY = 1;
    else
          KY = 1 - (LENY-1)*INCY;
              
              
    if (BETA!=1) 
        { if (INCY==1)
              {
              if (BETA==0) 
                  for(i= 0;i<LENY;i++)
                      Y[i] = 0;
              else
                  for(i = 1;i<LENY;i++)
                      Y[i] = BETA*Y[i];
              }  
              
          else
              { IY = KY;
                if(BETA==0) 
                      for(i=0;i<LENY;i++)
                        { Y[IY] = 0;
                          IY = IY + INCY;
                        } 
                else
                  for(i=0;i<LENY;i++)
                      { Y[IY]=BETA*Y[IY];
                        IY = IY + INCY;
                      }  
               }   
        }

              
    if (ALPHA==0) return ;
      if(TRANS=='N'||TRANS=='n') 
/*
*        Form  y := alpha*A*x + y.
*/
         { JX = KX;
           if(INCY==1)
              {for(j = 0;j<N;j++)
                 { if(X[JX]!=0) 
                      { TEMP = ALPHA*X[JX];
                         for(i=0;i<M;i++)
                          Y[i] = Y[i] + TEMP*A[i*LENX+j];
                      }
                   JX = JX + INCX;
                 }
              }   
          else
              { for(j =0;j<N;j++)
                  { if (X[JX]!=0) 
                      { TEMP = ALPHA*X[JX];
                        IY = KY;
                         for(i = 0;i<M;i++)
                          { Y[IY] = Y[IY] + TEMP*A[i*LENX+j];
                            IY = IY + INCY;
                          }
                       }      
                   JX = JX + INCX;
                  }
            }
          }
      else
      {
/*
*        Form  y := alpha*A'*x + y.
*/
          JY = KY;
           if(INCX==1)
             {for(j = 0;j<N;j++)
                  { TEMP = 0;
                    for(i = 0;i<M;i++)
                      TEMP = TEMP + A[i*LENX+j]*X[i];
                    Y[JY] = Y[JY] + ALPHA*TEMP;
                     JY = JY + INCY;
                  }   
             }         
          else
              {for(j=0;j<N;j++)
                  { TEMP = 0;
                    IX = KX;
                    for(i = 0;i<M;i++)
                       { TEMP = TEMP + A[i*LENX+j]*X[IX];
                         IX = IX + INCX;
                       }
                    Y[JY] = Y[JY] + ALPHA*TEMP;
                    JY = JY + INCY;
                  }
              }
                       
         } 
return;
}





