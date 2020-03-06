void STRMV(char UPLO,char TRANS,char DIAG,int N,float*A,int LDA,float* X,int INCX)
{
int KX,KY,IX,IY,JX,JY,i,j;
float TEMP;
/*
*     Quick return if possible.
*/
	if(N==0) return;
/*
      NOUNIT = LSAME(DIAG,'N')
*/
	if(INCX<0)
          KX = -(1)*(N-1)*INCX;
    else  if (INCX!=1) 
          KX = 0;
    if(TRANS=='N'||TRANS=='n') 
	  {	
/*
*        Form  x := A*x.
*/
	  	if(UPLO=='U'||UPLO=='u')
          {	if(INCX==1)
              {	for(j=0;j<N;j++)
                   {	if(X[j]!=0)
                          {	TEMP = X[j];
                           	for(i=0;i<=j-1;i++)
                              {	X[i] = X[i]+TEMP*A[i*N+j];
                              }
                            if(DIAG=='N'||DIAG=='n')
                             X[j] = X[j]*A[j*N+j];
						  }	                      
   		           }
   		      }     
            else
              {	JX = KX;
                for(j=0;j<N;j++)
                   {	if(X[JX]!=0)
                          {	TEMP = X[JX];
                          	IX = KX;
                          	for(i=0;i<=j-1;i++)
                              {	X[IX] = X[IX] + TEMP*A[i*N+j];
                              	IX = IX + INCX;
                              }	
                            if(DIAG=='N'||DIAG=='n')
                             X[JX] = X[JX]*A[j*N+j];
                          }
                      JX = JX + INCX;
                   }
              }
          } 
        else   
          {	if(INCX==1)
              {	for(j=N-1;j>=0;j--)
                   {	if(X[j]!=0)
                          {	TEMP = X[j];
                          	for(i=N-1;i>=j+1;i--)
                          	  {	X[i] = X[i] + TEMP*A[i*N+j];
                              }
                          	if(DIAG=='N'||DIAG=='n')
                          	 X[j] = X[j]*A[j*N+j];
 						  }	                     
  				   }
  			  }	   
            else
              {	KX = KX + (N-1)*INCX;
                JX = KX;
                for(j=N-1;j>=0;j--)
                   {	if(X[JX]!=0)
                          {	TEMP = X[JX];
                         	IX = KX;
                            for(i=N-1;i>=j+ 1;i--)
                              {	X[IX] = X[IX] + TEMP*A[i*N+j];
                              	IX = IX - INCX;
                              }	
                          	if(DIAG=='N'||DIAG=='n')
                          	  X[JX] = X[JX]*A[j*N+j];
                          }
                      JX = JX - INCX;
                   }
              }
          }
      }
      
    else
	  {
/*
*        Form  x := A'*x.
*/
       	if(UPLO=='U'||UPLO=='u')
          {	if(INCX==1)
              {	for(j=N-1;j>=0;j--)
                   {	TEMP = X[j];
                        if(DIAG=='N'||DIAG=='n')
                          TEMP = TEMP*A[j*N+j];
                        for(i =j-1;i>=0;i--)
                          TEMP = TEMP + A[i*N+j]*X[i];
                        X[j] = TEMP;
                   }
              }      
            else  
              {	JX = KX + (N-1)*INCX;
                for(j=N-1;j>=0;j--)
                   {	TEMP = X[JX];
                      	IX = JX;
                      	if(DIAG=='N'||DIAG=='n')
                      	  TEMP = TEMP*A[j*N+j];
                      	for(i =j-1;i>=0;i--)
                          {	IX = IX - INCX;
                          	TEMP = TEMP + A[i*N+j]*X[IX];
                          }  
                        X[JX] = TEMP;
                      	JX = JX - INCX;
                   }
              }
          }             
        else
          {	if(INCX==1)
              {	for(j=0;j<N;j++)
                   {	TEMP = X[j];
                        if(DIAG=='N'||DIAG=='n')
                          TEMP = TEMP*A[j*N+j];
                        for(i=j + 1;i<N;i++)
                          {	TEMP = TEMP + A[i*N+j]*X[i];
                          }
                      	X[j] = TEMP;
                   }
              }
            else
              {	JX = KX;
                for(j =0;j<N;j++)
                   {	TEMP = X[JX];
                      	IX = JX;
                      	if(DIAG=='N'||DIAG=='n')
                      	  TEMP = TEMP*A[j*N+j];
                      	for(i =j + 1;i<N;i++)
                          {	IX = IX + INCX;
                            TEMP = TEMP + A[i*N+j]*X[IX];
                          }
                        X[JX] = TEMP;
                        JX = JX + INCX;
                   }
              }
          }
      }

return;


}




