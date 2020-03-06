void SGBMV(char TRANS,int M,int N,int KL,int KU,float ALPHA,float* A,int LDA,float* X,int INCX,float BETA,float* Y,int INCY)
{

int LENX,LENY,KX,KY,KUP1,K,JX,JY,IX,IY,a,i,j;
float TEMP;



if((M==0)||(N==0)||((ALPHA==0)&&(BETA==1))) return ;


if(TRANS=='N'||TRANS=='n') 
  {	LENX = N;
    LENY = M;
  }
else
  {	LENX = M;
    LENY = N;
  }
if(INCX>0)
  KX = 0;
else
  KX =-1*(LENX-1)*INCX;
if(INCY>0)
  KY =0;
else
  KY = -1* (LENY-1)*INCY;


/*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through the band part of A.
*
*     First form  y := beta*y.
*/
if(BETA!=1)
   {	if(INCY==1)
          {if(BETA==0)
             {	for(i=0;i<LENY;i++)
                      Y[i] =0;
             }
           else
             {	for(i=0;i<LENY;i++)
                      Y[i] = BETA*Y[i];
             }
		  }
		else
          {	IY = KY;
            if(BETA==0) 
              {	for(i =0;i<LENY;i++)
                    {  Y[IY] = 0;
                       IY = IY + INCY;
                    }
              }
            else
              {	for(i =0;i<LENY;i++)
                    {  Y[IY] = BETA*Y[IY];
                       IY = IY + INCY;
                    }
              }
          }
   }   
if(ALPHA==0) return ;
     
KUP1 = KU + 1;
if(TRANS=='N'||TRANS=='n')
  {
/*
*        Form  y := alpha*A*x + y.
*/
   	JX = KX;
    if(INCY==1)
      {	for(j = 0;j<N;j++)
           { if(X[JX]!=0)
               {	TEMP = ALPHA*X[JX];
                    K = KUP1 - j;
                    a=(M<j+KL?M:j+KL);
                    for(i = (0>j-KU?0:j-KU);i<a;i++)
                          Y[i] = Y[i] + TEMP*A[(K+i)*N+j];//check
   	            }
   	         JX = JX + INCX;
          }
      }
    else
      {	for(j = 0;j<N;j++)
           { if(X[JX]!=0)
               {	TEMP = ALPHA*X[JX];
                    IY = KY;
                    K = KUP1 - j;
                    a=(M<j+KL?M:j+KL);
                    for(i =(0>j-KU?0:j-KU);i<a;i++)
                       {  Y[IY] = Y[IY] + TEMP*A[(K+i)*N+j];
                          IY = IY + INCY;
                       }   
               }
             JX = JX + INCX;
             if(j>KU) KY = KY + INCY;
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
      {	for(j =0;j<N;j++)
           { TEMP = 0;
             K = KUP1 - j;
             a=(M<j+KL?M:j+KL);
             for(i=(0>j-KU?0:j-KU);i<a;i++)
                  TEMP = TEMP + A[(K+i)*N+j]*X[i];
             Y[JY] = Y[JY] + ALPHA*TEMP;
             JY = JY + INCY;
           }
      }     
    else
      {for(j =0;j<N;j++)
          { TEMP =0;
            IX = KX;
            K = KUP1 - j;
            a=(M<j+KL?M:j+KL);
            for(i=(0>j-KU?0:j-KU);i<a;i++)
               { TEMP = TEMP + A[(K+i)*N+j]*X[IX];
                 IX = IX + INCX;
               }
            Y[JY] = Y[JY] + ALPHA*TEMP;
            JY = JY + INCY;
            if(j>KU) KX = KX + INCX;
          } 
     }
  } 
return;  

}
