void SSYMV(char UPLO,int N,float ALPHA,float* A,int LDA,float* X,int INCX,float BETA,float* Y,int INCY)
{
	int i,j,KX,KY,JX,JY,IX,IY;
    float TEMP1,TEMP2;

/* Quick Return */

	      if (((N==0)||(ALPHA==0))&&(BETA==1)) return ;
	      
	      if(INCX>0) 
              KX = 0;
      	  else
              KX = -1* (N-1)*INCX;
          if (INCY>0) 
              KY = 0;
          else
              KY = -1*(N-1)*INCY;
              
              
/*                    Calculation of Y = Y*BETA                        */      
     
      if(BETA!=1) 
          {	if(INCY==1)
              {	if(BETA==0)
                  {	for(i = 0;i<N;i++)
                       Y[i]= 0;
                  }
                else
                  {	for(i = 0;i<N;i++)
                       Y[i] = BETA*Y[i];
                  }
              }
            else
              {	IY = KY;
              	if(BETA==0)
              	  {	for(i = 0;i<N;i++)
                      {	Y[IY] = 0;
                      	IY = IY + INCY;
                      }
                  }              
                else
                  {	for(i = 0;i<N;i++)
                      {	Y[IY] = BETA*Y[IY];
                       	IY = IY + INCY;
                      }
                  }    
              }
          } 
           
/* when ALPHA is zero output is Y*BETA which has already been calcuated in Y.*/
           
if(ALPHA==0) return ;   


if(UPLO=='U'||UPLO=='u') 
{
/*
       Form  y  when A is stored in upper triangle.
*/
	if((INCX==1)&&(INCY==1))
      {	for(j = 0;j<N;j++)
           {	TEMP1 = ALPHA*X[j];
                TEMP2 = 0;
                for(i = 0;i<=j-1;i++)
                  {   Y[i] = Y[i] + TEMP1*A[i*N+j];
                      TEMP2 = TEMP2 + A[i*N+j]*X[i];
   		          }
                  Y[j] = Y[j] + TEMP1*A[j*N+j] + ALPHA*TEMP2;
           }       
      }
    else
      {	JX = KX;
        JY = KY;
        for(j = 0;j<N;j++)
           {	TEMP1 = ALPHA*X[JX];
                TEMP2 = 0;
                IX = KX;
                IY = KY;
                for(i = 0;i<=j - 1;i++)
                   {	Y[IY] = Y[IY] + TEMP1*A[i*N+j];
                      	TEMP2 = TEMP2 + A[i*N+j]*X[IX];
                        IX = IX + INCX;
                        IY = IY + INCY;
                   }    
                Y[JY] = Y[JY] + TEMP1*A[j*N+j] + ALPHA*TEMP2;
                JX = JX + INCX;
                JY = JY + INCY;
           }
      }
}          
else
{
/*
*        Form  y  when A is stored in lower triangle.
*/
    if((INCX==1)&&(INCY==1))
      {	for(j = 0;j<N;j++)
           {	TEMP1 = ALPHA*X[j];
                TEMP2 = 0;
                Y[j] = Y[j] + TEMP1*A[j*N+j];
                for(i = j+1;i<N;i++)
                   {	Y[i] = Y[i] + TEMP1*A[i*N+j];
                      	TEMP2 = TEMP2 + A[i*N+j]*X[i];
                   }
                Y[j] = Y[j] + ALPHA*TEMP2;
  	       }
  	  }     	
     else
      {	JX = KX;
        JY = KY;
        for(j=0;j<N;j++)
           {	TEMP1 = ALPHA*X[JX];
                  TEMP2 = 0;
                  Y[JY] = Y[JY] + TEMP1*A[j*N+j];
                  IX = JX;
                  IY = JY;
                  for(i = j + 1;i<N;i++)
                     {	IX = IX + INCX;
                      	IY = IY + INCY;
                      	Y[IY] = Y[IY] + TEMP1*A[i*N+j];
                      	TEMP2 = TEMP2 + A[i*N+j]*X[IX];
                     } 
                  Y[JY] = Y[JY] + ALPHA*TEMP2;
                  JX = JX + INCX;
                  JY = JY + INCY;
           }
      }
}
return ;
}  
