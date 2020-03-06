#include "blas2.h"
#include "tbb/tbb.h"
#include <cstdio>
#include <iostream>

using namespace std;
using namespace tbb;

class MatrixMultiplyBody2D
{
    float *my_A, *my_X, *my_Y;
    float my_ALPHA;
    int my_ROWS,my_COLS,my_INCX,my_INCY;

public:
    void operator()( const blocked_range2d<size_t>& r ) const
    {
#if SHOW_THREAD_INFO == 1
        cout<<"\n TRNS range: "<<r.rows().begin()<<","<<r.cols().begin()<<" to "<<r.rows().end()<<","<<r.cols().end()<<" object: "<<this<<", thead-id: "<<this_tbb_thread::get_id();
#endif
        int cols = my_COLS;
        int rows = my_ROWS;
        float *a = my_A;



        float *x = my_X;
        int incX=my_INCX;
        float alpha=my_ALPHA;

        int firstIndexX=0;
        if (incX < 0)
            firstIndexX=(cols-1)*incX*-1;
        int indexX=firstIndexX+r.rows().begin()*incX;


        float *y = my_Y;
        int incY=my_INCY;

        int firstIndexY=0;
        if (incY < 0)
            firstIndexY=(rows-1)*incY*-1;
        firstIndexY=firstIndexY+r.cols().begin()*incY;
        int indexY=firstIndexY;

        int indexA=r.rows().begin()*cols + r.cols().begin();


        for ( size_t i=r.rows().begin(); i!=r.rows().end(); ++i,indexX+=incX )
        {
            for ( size_t j=r.cols().begin(); j!=r.cols().end(); ++j, indexY+=incY )
            {
//              if(indexA==250)
//                cout<<alpha<<"* x["<<indexX<<"]("<<x[indexX]<<")*y["<<indexY<<"]("<<y[indexY]<<")+ a["<<indexA<<"]("<<a[indexA]<<")";
                a[indexA]+= alpha*x[indexX]*y[indexY];
                indexA++;
            }
            indexA+=cols-r.cols().end()+r.cols().begin();
            indexY=firstIndexY;
        }
    }
    MatrixMultiplyBody2D( int rows,int cols,float *a, float *x, float *y,float alpha, int incX, int incY ) :
            my_ROWS(rows), my_COLS(cols), my_A(a), my_X(x), my_Y(y), my_ALPHA(alpha), my_INCX(incX), my_INCY(incY) {}
};


//     A := alpha*x*y' + A
void SGER(int M,int N,float ALPHA,float* X,int INCX,float* Y,int INCY,float* A,int LDA)
{
    parallel_for (/*Range*/ blocked_range2d<size_t>(0, M, TILESIZE, 0,N, TILESIZE), MatrixMultiplyBody2D(M,N,A,X,Y,ALPHA,INCX,INCY) );

    return;
}

