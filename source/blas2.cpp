#include <iostream>
#include <cstdio>
#include <iomanip>
#include <ctime>
#include <cstdlib>
#include "blas2.h"
#include "tbb/tick_count.h"


using namespace std;
using namespace tbb;

void generate_GBM(float *A,int KU,int KL,int rows,int cols,float *x,float *y,float *res,float *res1,float &alpha, float &beta);
void intializeVars_rank(float *A,float *x,float *y,float &alpha,int rows, int cols);
void intializeVars(float *,float *,float *,float *, float*,float &, float &, int, int);
void intializeVars(float *,float *,float *, int);
void printArrays_rank(float *A,float *x,float *y,float *res,int, int);
void printArrays(float *,float *,float *,float *,int, int);
void printArrays(float *A,float *x,float *res,int rows);
void printArray(float *,char const *,int, int);
void copyArray(float *,float *,int, int);
float* transpose(float *,int,int);
bool checkForEquality(float *,float *,int, int);

void performTriangleMultiplication(int rows,int incX,int tileSize,int verify,int echo);
float* normal_function(char upper_lower,char diag_zero,int rows,float *A,float *x,int incX,float *res);
void normal_function(int ,int ,float ,float* ,int ,float* ,int ,float ,float*,int );
void normal_function_rank(int rows,int cols,float alpha,float *x,int incX,float *y,int incY,float *res1);

void ssymv(char t,int rows,int cols,float alpha,float *A,int lda,float *x,int incX,float beta,float *y,int incY)
{
    SSYMV(t=='n'?'u':'l',rows,alpha,A,rows,x,incX,beta,y,incY);
}

void strmv(char t,int rows,int cols,float alpha,float *A,int lda,float *x,int incX,float beta,float *y,int incY)
{
    //  STRMV(uplo,'n',rows,alpha,A,rows,x,incX);
}

void sger(char t,int rows,int cols,float alpha,float *A,int lda,float *x,int incX,float beta,float *y,int incY)
{
    //  STRMV(uplo,'n',rows,alpha,A,rows,x,incX);
}

void sgbmv(char t,int rows,int cols,float alpha,float *A,int lda,float *x,int incX,float beta,float *y,int incY)
{
    //  STRMV(uplo,'n',rows,alpha,A,rows,x,incX);
}

void strmvTD(char uplo,char t,char d,int rows,float *A,int lda,float *x,int incX)
{
    STRMV(uplo,t,d,rows,A,rows,x,incX);
}

int lenX,lenY=0;
int incX=DEF_INCX;
int incY=DEF_INCY;

string func_names[]={"SGEMV","SSYMV","STRMV","SGBMV", "SGER"};
int func_count=5;
string func_name=func_names[0];
int func_id=0;

int main(int argc,char *argv[])
{

    void (*blas_function)(char,int,int,float,float* ,int ,float* ,int ,float ,float* ,int );
    int rows=DEF_ROWS;
    int cols=DEF_COLS;
    int tileSize=NO_TILES;
    int verify=DEF_VERIFY;
    int echo=DEF_ECHO;

    void (*blas_funcns[]) (char,int,int,float,float* ,int ,float* ,int ,float ,float* ,int )={SGEMV,ssymv,strmv,sgbmv,sger};
    blas_function=blas_funcns[0];

    cout<<"Argument list: function_name rows cols incX incY tileSize verify(1 or 0) echo(1 or 0)"<<endl;

    if (argc>=2)
    {
        int i=0;
        for (i=0;i<func_count;i++)
            if (func_names[i]==argv[1])
            {
                blas_function=blas_funcns[i];
                func_name=func_names[i];
                func_id=i;
                break;
            }

        if (i==func_count)
        {
            cout<<"Wrong choice of blas function. Instead of "<<argv[1]<<" use";
            for (int i=0;i<func_count;i++)
                cout<<" "<<func_names[i];
            cout<<endl;
            return -1;
        }

        if (argc>=3)
        {
            sscanf(argv[2],"%d",&rows);
            if (argc>=4)
            {
                sscanf(argv[3],"%d",&cols);
                if (argc>=5)
                {
                    sscanf(argv[4],"%d",&incX);
                    if (argc>=6)
                    {
                        sscanf(argv[5],"%d",&incY);
                        if (argc>=7)
                        {
                            sscanf(argv[6],"%d",&tileSize);
                            if (argc>=8)
                            {
                                sscanf(argv[7],"%d",&verify);
                                if (argc==9)
                                    sscanf(argv[8],"%d",&echo);
                            }
                        }
                    }
                }
            }
        }
        if (func_id==1 || func_id==2)
        {
            cols=rows;	//ssymv or strmv
        }

    }

    srand(100);
    switch (func_id)
    {
    case 4:
    {
        float *A=NULL,*x=NULL,*y=NULL,*res=NULL,*res1=NULL;
        float alpha=0;
        try
        {
            A=new float[rows*cols];

            lenX=incX<0?1+(cols-1)*incX*-1:1+(cols-1)*incX;
            x=new float[lenX];

            lenY=incY<0?1+(rows-1)*incY*-1:1+(rows-1)*incY;
            y=new float[lenY];

            intializeVars_rank(A,x,y,alpha,rows,cols);

            res=new float[rows*cols];
            copyArray(A,res,rows,cols);

            tick_count start,end;

            start =tick_count::now();
#if NORMAL == 1
            normal_function_rank(rows,cols,alpha,x,incX,y,incY,res);
#else
            SGER(rows,cols,alpha,x,incX,y,incY,res,0);
#endif
            end=tick_count::now();

            cout<<func_name<<" , alpha:"<<alpha<<" rows:"<<rows<<" , cols:"<<cols<<" ,time:"<<(end-start).seconds()<<endl;

            if (echo)
                printArrays_rank(A,x,y,res,rows,cols);

            if (verify)
            {
                res1=new float[rows*cols];
                copyArray(A,res1,rows,cols);
                normal_function_rank(rows,cols,alpha,x,incX,y,incY,res1);
                if (!checkForEquality(res,res1,rows,cols))
                {
                    cout<<endl<<"problem with "<<func_name<<endl;
                    if (echo)	printArray(res1,"normal res",rows,cols);
                    throw -1;
                }

                delete []res1;
                delete []res;
                res1=res=NULL;
            }
        }
        catch (...)
        {
            if (echo)
                cout<<endl<<"freeing memory"<<endl;

            if (A!=NULL) delete []A;
            if (x!=NULL) delete []x;
            if (y!=NULL) delete []y;
            if (res1!=NULL) delete []res1;
            if (res!=NULL) delete []res;
            A=x=y=res1=res=NULL;
            throw;
        }

        if (echo)
            cout<<endl<<"freeing memory"<<endl;

        if (A!=NULL) delete []A;
        if (x!=NULL) delete []x;
        if (y!=NULL) delete []y;
        if (res1!=NULL) delete []res1;
        if (res!=NULL) delete []res;
        A=x=y=res1=res=NULL;


    }
    break;
    case 3:
    {
        float *A,*x,*y,*res,*res1,*res2=NULL;
        int KU,KL;
        KU=rand()%rows;
        KL=rand()%cols;
        A=new float[rows*cols];

        lenX=incX<0?1+(cols-1)*incX*-1:1+(cols-1)*incX;
        x=new float[lenX];

        lenY=incY<0?1+(rows-1)*incY*-1:1+(rows-1)*incY;
        y=new float[lenY];

        res=new float[lenY];
        res1=new float[lenX];
        float alpha,beta=0;
        alpha=0;
        beta=0;
        generate_GBM(A,KU,KL,rows,cols,x,y,res,res1,alpha,beta);
        tick_count start,end;

        if (echo) printArrays(A,x,res,rows);
        start =tick_count::now();
        SGBMV('n',rows,cols,KL,KU,alpha,A,rows,x,incX,beta,y,incY);
        cout<<func_name<<" , alpha:"<<alpha<<" beta:"<<beta<<" rows:"<<rows<<" , cols:"<<cols<<",KU = "<<KU<<",KL = "<<KL<< ",time:"<<(end-start).seconds()<<endl;
        if (verify)
        {
            res2=new float[lenY];
            copyArray(y,res2,1,lenY);
            normal_function(rows,cols,alpha,A,rows,x,incX,beta,res2,incY);
            if (!checkForEquality(res,res2,1,lenY))
            {
                cout<<endl<<"problem with "<<func_name<<endl;
                if (echo)	printArray(res2,"normal res",1,lenY);
                throw -1;
            }
        }
    }
    break;
    case 2:
        performTriangleMultiplication(rows,incX,tileSize,verify,echo);
        break;
    default:
    {
        float *A=NULL,*x=NULL,*y=NULL,*res=NULL,*res1=NULL,*res2=NULL;

        try
        {
            A=new float[rows*cols];

            lenX=incX<0?1+(cols-1)*incX*-1:1+(cols-1)*incX;
            x=new float[lenX];

            lenY=incY<0?1+(rows-1)*incY*-1:1+(rows-1)*incY;
            y=new float[lenY];

            res=new float[lenY];
            res1=new float[lenX];
            float alpha,beta=0;
            alpha=0;
            beta=0;

            intializeVars(A,x,y,res,res1,alpha,beta,rows,cols);

            tick_count start,end;

            start =tick_count::now();
#if NORMAL == 1
            normal_function(rows,cols,alpha,A,rows,x,incX,beta,res,incY);
#else
            blas_function('n',rows,cols,alpha,A,rows,x,incX,beta,res,incY);
#endif
            end=tick_count::now();

            cout<<func_name<<" , alpha:"<<alpha<<" beta:"<<beta<<" rows:"<<rows<<" , cols:"<<cols<<" ,time:"<<(end-start).seconds()<<endl;

            if (echo)
                printArrays(A,x,y,res,rows,cols);

            if (verify)
            {
                res2=new float[lenY];
                copyArray(y,res2,1,lenY);
                normal_function(rows,cols,alpha,A,rows,x,incX,beta,res2,incY);
                if (!checkForEquality(res,res2,1,lenY))
                {
                    cout<<endl<<"problem with "<<func_name<<endl;
                    if (echo)	printArray(res2,"normal res",1,lenY);
                    throw -1;
                }

                delete []res2;
                delete []res;
                res2=res=NULL;
            }
#if TRANSPOSE == 1
            start =tick_count::now();
            blas_function('t',rows,cols,alpha,A,rows,y,incY,beta,res1,incX);
            end=tick_count::now();

            cout<<func_name<<" , rows:"<<rows<<" , cols:"<<cols<<" ,transposed, time:"<<(end-start).seconds()<<endl;
            if (echo)	printArray(res1,"blas Trans_res",1,lenX);

            if (verify)
            {
                res2=new float[lenX];
                copyArray(x,res2,1,lenX);
                A=transpose(A,rows,cols);
                if (echo) printArray(A,"trans_a",cols,rows);
                normal_function(cols,rows,alpha,A,cols,y,incY,beta,res2,incX);
                if (!checkForEquality(res1,res2,1,lenX))
                {
                    cout<<endl<<"problem in trans with "<<func_name<<endl;
                    if (echo)	printArray(res2,"normal trans_res",1,lenX);
                    throw -1;
                }
            }
#endif
        }
        catch (...)
        {
            if (echo)
                cout<<endl<<"freeing memory"<<endl;

            if (A!=NULL) delete []A;
            if (x!=NULL) delete []x;
            if (y!=NULL) delete []y;
            if (res1!=NULL) delete []res1;
            if (res2!=NULL) delete []res2;
            if (res!=NULL) delete []res;
            A=x=y=res1=res2=res=NULL;
            throw;
        }

        if (echo)
            cout<<endl<<"freeing memory"<<endl;

        if (A!=NULL) delete []A;
        if (x!=NULL) delete []x;
        if (y!=NULL) delete []y;
        if (res1!=NULL) delete []res1;
        if (res2!=NULL) delete []res2;
        if (res!=NULL) delete []res;
        A=x=y=res1=res2=res=NULL;
    }
    }
    return 0;
}

void performTriangleMultiplication(int rows,int incX,int tileSize,int verify,int echo)
{
    float *A=NULL,*x=NULL,*res=NULL,*res1=NULL,*res2=NULL;

    try
    {
        A=new float[rows*rows];

        lenX=incX<0?1+(rows-1)*incX*-1:1+(rows-1)*incX;
        x=new float[lenX];

        res=new float[lenX];

        intializeVars(A,x,res,rows);

        tick_count start,end;
        string func_name="strmv";

        /******* x=u(A)*x (Diagonal elements included) *******/
        start =tick_count::now();
#if NORMAL == 1
        normal_function('u','n',rows,A,x,incX,res1);
#else
        strmvTD('u','n','n',rows,A,rows,res,incX);
#endif
        end=tick_count::now();
        cout<<func_name<<" ,rows:"<<rows<<" ,upper non-transposed time:"<<(end-start).seconds()<<endl;

        if (echo) printArrays(A,x,res,rows);

        if (verify)
        {
            res1=new float[lenX];
            copyArray(x,res1,1,lenX);
            normal_function('u','n',rows,A,x,incX,res1);
            if (!checkForEquality(res,res1,1,lenX))
            {
                cout<<endl<<"problem with "<<func_name<<endl;
                if (echo)	printArray(res1,"normal res",1,lenX);
                throw -1;
            }

        }


//        /******* x=l(A')*x (Diagonal elements included) *******/
//        copyArray(x,res,1,lenX);
//        start =tick_count::now();
//        strmvTD('u','t','n',rows,A,rows,res,incX);
//        end=tick_count::now();
//        cout<<endl<<"lower transposed time:"<<(end-start).seconds()<<endl;
//
//        transpose(A,rows,rows);
//        if (echo) printArrays(A,x,res,rows);
//
//        if (verify)
//        {
//            res2=new float[lenX];
//            copyArray(x,res2,1,lenX);
//            normal_function('l','n',rows,A,x,incX,res2);
//            if (!checkForEquality(res,res2,1,lenX))
//            {
//                cout<<endl<<"problem with "<<func_name<<endl;
//                if (echo)	printArray(res2,"normal res",1,lenX);
//                throw -1;
//            }
//
//        }
//
//        /******* x=l(A)*x (Diagonal elements included) *******/
//        copyArray(x,res,1,lenX);
//        start =tick_count::now();
//        strmvTD('l','n','n',rows,A,rows,res,incX);
//        end=tick_count::now();
//        cout<<endl<<"lower non-transposed time:"<<(end-start).seconds()<<endl;
//
//        if (echo) printArray(res,"blas res",1,rows);
//
//        if (verify)
//        {
//            if (!checkForEquality(res,res2,1,lenX))
//            {
//                cout<<endl<<"problem with "<<func_name<<endl;
//                if (echo)	printArray(res2,"normal res",1,lenX);
//                throw -1;
//            }
//
//        }
//
//
//        /******* x=u(A')*x (Diagonal elements included) *******/
//        copyArray(x,res,1,lenX);
//        start =tick_count::now();
//        strmvTD('l','t','n',rows,A,rows,res,incX);
//        end=tick_count::now();
//        cout<<endl<<"upper transposed time:"<<(end-start).seconds()<<endl;
//
//        if (echo) printArray(res,"blas res",1,rows);
//
//        if (verify)
//        {
//            if (!checkForEquality(res,res1,1,lenX))
//            {
//                cout<<endl<<"problem with "<<func_name<<endl;
//                if (echo)	printArray(res1,"normal res",1,lenX);
//                throw -1;
//            }
//
//        }
//
//        //Make all diagonal elements 0
//        for (int i=0;i<rows;i++)
//            A[rows*i + i]=1;
//
//
//        /******* x=l(A)*x (Diagonal elements 0) *******/
//        copyArray(x,res,1,lenX);
//        start =tick_count::now();
//        strmvTD('l','n','y',rows,A,rows,res,incX);
//        end=tick_count::now();
//        cout<<endl<<"******** All diagonal elements made 0 ********"<<endl<<"lower non-transposed time:"<<(end-start).seconds()<<endl;
//
//        if (echo) printArrays(A,x,res,rows);
//
//        if (verify)
//        {
//            copyArray(x,res1,1,lenX);
//            normal_function('l','y',rows,A,x,incX,res1);
//            if (!checkForEquality(res,res1,1,lenX))
//            {
//                cout<<endl<<"problem with "<<func_name<<endl;
//                if (echo)	printArray(res1,"normal res",1,lenX);
//                throw -1;
//            }
//
//        }
//
//
//        /******* x=u(A')*x (Diagonal elements 0) *******/
//        copyArray(x,res,1,lenX);
//        start =tick_count::now();
//        strmvTD('l','t','y',rows,A,rows,res,incX);
//        end=tick_count::now();
//        cout<<endl<<"upper transposed time:"<<(end-start).seconds()<<endl;
//
//        transpose(A,rows,rows);
//        if (echo) printArrays(A,x,res,rows);
//
//        if (verify)
//        {
//            copyArray(x,res2,1,lenX);
//            normal_function('u','y',rows,A,x,incX,res2);
//            if (!checkForEquality(res,res2,1,lenX))
//            {
//                cout<<endl<<"problem with "<<func_name<<endl;
//                if (echo)	printArray(res2,"normal res",1,lenX);
//                throw -1;
//            }
//
//        }
//
//        /******* x=u(A)*x (Diagonal elements 0) *******/
//        copyArray(x,res,1,lenX);
//        start =tick_count::now();
//        strmvTD('u','n','y',rows,A,rows,res,incX);
//        end=tick_count::now();
//        cout<<endl<<"upper non-transposed time:"<<(end-start).seconds()<<endl;
//
//        if (echo) printArray(res,"blas res",1,rows);
//
//        if (verify)
//        {
//            if (!checkForEquality(res,res2,1,lenX))
//            {
//                cout<<endl<<"problem with "<<func_name<<endl;
//                if (echo)	printArray(res2,"normal res",1,lenX);
//                throw -1;
//            }
//
//        }
//
//
//        /******* x=l(A')*x (Diagonal elements 0) *******/
//        copyArray(x,res,1,lenX);
//        start =tick_count::now();
//        strmvTD('u','t','y',rows,A,rows,res,incX);
//        end=tick_count::now();
//        cout<<endl<<"lower transposed time:"<<(end-start).seconds()<<endl;
//
//        if (echo) printArray(res,"blas res",1,rows);
//
//        if (verify)
//        {
//            if (!checkForEquality(res,res1,1,lenX))
//            {
//                cout<<endl<<"problem with "<<func_name<<endl;
//                if (echo)	printArray(res1,"normal res",1,lenX);
//                throw -1;
//            }
//
//        }

    }
    catch (...)
    {
        if (echo)
            cout<<endl<<"freeing memory"<<endl;
        if (A!=NULL) delete []A;
        if (x!=NULL) delete []x;
        if (res1!=NULL) delete []res1;
        if (res2!=NULL) delete []res2;
        if (res!=NULL) delete []res;
        A=x=res1=res2=res=NULL;
        throw;
    }

    if (echo)
        cout<<endl<<"freeing memory"<<endl;
    if (A!=NULL) delete []A;
    if (x!=NULL) delete []x;
    if (res1!=NULL) delete []res1;
    if (res2!=NULL) delete []res2;
    if (res!=NULL) delete []res;
    A=x=res1=res2=res=NULL;
}


void intializeVars(float *A,float *x,float *res, int rows)
{
    float minimum=2.0;
    if (ALLOW_ZEROES)
        minimum=0.0;

    init_trng(A,rows,minimum);

    for (int i=0;i<lenX;i++)
    {
        x[i]=rand()%VALUE_RANGE + minimum;
        res[i]=x[i];
    }

}


void intializeVars_rank(float *A,float *x,float *y,float &alpha,int rows, int cols)
{
    float minimum=2.0;
    if (ALLOW_ZEROES)
        minimum=0.0;
    init_gen(A,rows,cols,minimum);

    for (int i=0;i<lenX;i++)
    {
        x[i]=rand()%VALUE_RANGE + minimum;
    }

    for (int i=0;i<lenY;i++)
    {
        y[i]=rand()%VALUE_RANGE + minimum;
    }

    alpha=rand()%VALUE_RANGE + minimum;

}

void intializeVars(float *A,float *x,float *y,float *res,float *res1,float &alpha, float &beta,int rows, int cols)
{
    float minimum=2.0;
    if (ALLOW_ZEROES)
        minimum=0.0;
    switch (func_id)
    {

    case 0:
        init_gen(A,rows,cols,minimum);
        break;
    case 1:
        init_sym(A,rows,minimum);
        break;
    }

    for (int i=0;i<lenX;i++)
    {
        x[i]=rand()%VALUE_RANGE + minimum;
        res1[i]=x[i];
    }

    for (int i=0;i<lenY;i++)
    {
        y[i]=rand()%VALUE_RANGE + minimum;
        res[i]=y[i];
    }

    alpha=rand()%VALUE_RANGE + minimum;
    beta=rand()%VALUE_RANGE + minimum;

}

void printArrays_rank(float *A,float *x,float *y,float *res,int rows, int cols)
{
    printArray(A,"A",rows,cols);
    printArray(x,"x",1,lenX);
    printArray(y,"y",1,lenY);
    printArray(res,"blas res",rows,cols);
}

void printArrays(float *A,float *x,float *y,float *res,int rows,int cols)
{
    printArray(A,"A",rows,cols);
    printArray(x,"x",1,lenX);
    printArray(y,"y",1,lenY);
    printArray(res,"blas res",1,rows);
}

void printArrays(float *A,float *x,float *res,int rows)
{
    printArray(A,"A",rows,rows);
    printArray(x,"x",1,lenX);
    printArray(res,"blas res",1,rows);
}

void printArray(float *Ar,char const *c,int rows,int cols)
{
    cout<<endl<<c<<"=";
    for (int i=0;i<rows;i++)
    {
        cout<<endl;
        for (int j=0;j<cols;j++)
        {
            cout<<" "<<Ar[i*cols + j];
        }
    }
    cout<<endl;
}

void copyArray(float *from,float *to,int rows,int cols)
{
    for (int i=0;i<rows;i++)
    {
        for (int j=0;j<cols;j++)
        {
            to[i*cols + j]=from[i*cols + j];
        }
    }
}

bool checkForEquality(float *a1,float *a2,int rows, int cols)
{
    for (int i=0;i<rows;i++)
    {
        for (int j=0;j<cols;j++)
        {
            if (a1[i*cols + j]!=a2[i*cols + j])
            {
                cout<<"\nValues differ at ("<<i<<","<<j<<") - parallel:"<<a1[i*cols + j]<<"  normal:"<<a2[i*cols + j]<<"\n";
                return false;
            }
        }
    }
    return true;
}

float* transpose(float *A,int rows,int cols)
{
    if (rows == cols)
    {
        for (int i=0;i<rows;i++)
            for (int j=i+1;j<rows;j++)
            {
                int temp=A[rows*i + j];
                A[rows*i + j]=A[rows*j + i];
                A[rows*j + i]=temp;
            }

        return A;
    }
    else
    {
        float *b=new float[rows*cols];
        int index=0;

        for (int i=0;i<cols;i++)
            for (int j=i;j<cols*rows;j+=cols,index++)
                b[index]=A[j];

        delete []A;
        return b;
    }


}

float* normal_function(char upper_lower,char diag_zero,int rows,float *A,float *x,int incX,float *res)
{

    int indexX=0;
    if (incX < 0)
        indexX=(rows-1)*incX*-1;
    int temp=indexX;
    int resIndex=indexX;

    if (diag_zero=='n')
    {
        int indexA=0;
        if (upper_lower=='l'||upper_lower=='L')
        {

            for (int i=0;i<rows;i++,resIndex+=incX,indexX=temp)
            {
                float value=0;

                for (int j=0;j<=i;j++,indexX+=incX)
                    value+=A[indexA++] * x[indexX];

                indexA=indexA+rows-i-1;
                res[resIndex]=value;
            }
        }
        else
        {

            for (int i=0;i<rows;i++,resIndex+=incX,indexX=temp)
            {
                float value=0;

                for (int j=i;j<rows;j++,indexX+=incX)
                    value+=A[indexA++] * x[indexX];

                indexA=indexA+i+1;
                res[resIndex]=value;
                temp+=incX;
            }
        }
    }
    else
    {
        if (upper_lower=='l'||upper_lower=='L')
        {
            int indexA=rows;

            for (int i=1;i<rows;i++,resIndex+=incX,indexX=temp)
            {
                float value=0;

                for (int j=0;j<i;j++,indexX+=incX)
                    value+=A[indexA++] * x[indexX];

                indexA=indexA+rows-i;
                res[resIndex]=value;
            }
        }
        else
        {
            int indexA=1;
            for (int i=0;i<rows;i++,resIndex+=incX,indexX=temp)
            {
                float value=0;

                for (int j=i+1;j<rows;j++,indexX+=incX)
                    value+=A[indexA++] * x[indexX];

                indexA=indexA+i+2;
                res[resIndex]=value;
            }
        }
    }

    return res;
}

void normal_function(int rows,int cols,float alpha,float *A,int lda,float *x,int incX,float beta,float *y,int incY)
{

    int indexY=0;
    if (incY < 0)
        indexY=(rows-1)*incY*-1;

    int indexA=0;
    for (int i=0;i<rows;i++,indexY+=incY)
    {
        float value=0;

        int indexX=0;
        if (incX < 0)
            indexX=(cols-1)*incX*-1;

        for (int j=0;j<cols;j++,indexX+=incX)
            value+=A[indexA++] * x[indexX];

        y[indexY]=alpha*value + beta*y[indexY];
    }

    if (indexA!=rows*cols)
    {
        cout<<"********* Error in multiplication. Access of values outside array A. ********";
    }


}

void normal_function_rank(int rows,int cols,float alpha,float *x,int incX,float *y,int incY,float *A)
{

    int indexX=0;
    if (incX < 0)
        indexX=(cols-1)*incX*-1;

    int firstIndexY=0;
    if (incY < 0)
        firstIndexY=(rows-1)*incY*-1;
    int indexY=firstIndexY;

    int indexA=0;
    for (int i=0;i<rows;i++,indexX+=incX)
    {
        for (int j=0;j<cols;j++,indexY+=incY)
        {
//            if(indexA==250)
//              cout<<alpha<<"* x["<<indexX<<"]("<<x[indexX]<<")*y["<<indexY<<"]("<<y[indexY]<<")+ a["<<indexA<<"]("<<A[indexA]<<")";
            A[indexA]= y[indexY]*alpha*x[indexX] + A[indexA];
            indexA++;
        }
        indexY=firstIndexY;
    }

    if (indexA!=rows*cols)
    {
        cout<<"********* Error in multiplication. Access of values outside array A. ********";
    }

}

void generate_GBM(float *A,int KU,int KL,int rows,int cols,float *x,float *y,float *res,float *res1,float &alpha, float &beta)
{
    int i,j;
    float minimum=2.0;
    if (ALLOW_ZEROES)
        minimum=0.0;
    for (i=0;i<rows;i++)
    {
        for (j=0;j<cols;j++)
        {
            if (j<KU+i&&i<KL+j)
                A[i*rows + j]=rand()%VALUE_RANGE + minimum;
            else
                A[i*rows+j]=0;
        }
    }
    for (int i=0;i<lenX;i++)
    {
        x[i]=rand()%VALUE_RANGE + minimum;
        res1[i]=x[i];
    }

    for (int i=0;i<lenY;i++)
    {
        y[i]=rand()%VALUE_RANGE + minimum;
        res[i]=y[i];
    }

    alpha=rand()%VALUE_RANGE + minimum;
    beta=rand()%VALUE_RANGE + minimum;

}
