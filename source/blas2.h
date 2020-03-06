#define VALUE_RANGE 50
#define NO_TILES 0
#define TILESIZE 256
#define DEF_COLS 500
#define DEF_ROWS 500
#define DEF_VERIFY 0
#define DEF_ECHO 0
#define DEF_INCX 1
#define DEF_INCY 1
#define ALLOW_ZEROES false

#define min(X, Y)  ((X) < (Y) ? (X) : (Y))

#define SHOW_THREAD_INFO 0
#define SHOW_INNER_THREAD_INFO 0
#define INNER_TIME 0

#define NORMAL 1
#define TRANSPOSE 0

#define PARTITION 0   /*0: auto, 1: simple, 2:affinity*/
#define GRAINSIZE 1000

#define REDUCE 0      /*1: use reduction, 0: otherwise */

#define USE_TEMP_ARR 0   /*0: use temporary array, 1: otherwise */

#define USE_TASKS 1
#define CUTOFF 25000
#define PROCESSORS 4



void SGEMV(char TRANS,int M,int N,float ALPHA,float* A,int LDA,float* X,int INCX,float BETA,float* Y,int INCY);
void SSYMV(char UPLO,int N,float ALPHA,float* A,int LDA,float* X,int INCX,float BETA,float* Y,int INCY);
void STRMV(char UPLO,char TRANS,char DIAG,int N,float*A,int LDA,float* X,int INCX);
void SGBMV(char TRANS,int M,int N,int KL,int KU,float ALPHA,float* A,int LDA,float* X,int INCX,float BETA,float* Y,int INCY);
void SGER(int M,int N,float ALPHA,float* X,int INCX,float* Y,int INCY,float* A,int LDA);
void init_gen(float *A, int rows, int cols, float minimum);
void init_sym(float *A, int rows, float minimum);
void init_trng(float *A, int rows, float minimum);
