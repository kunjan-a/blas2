#include "blas2.h"
#include "tbb/tbb.h"
#include <cstdio>
#include <iostream>

using namespace std;
using namespace tbb;


class ChunkTask_GEN : public task
//class ChunkTask_GEN
{
public:
    int my_firstRow,my_CHUNKSIZE,my_ROWS,my_COLS,my_KX,my_KY,my_INCX,my_INCY;
    float *my_A,*my_X,*my_Y;
    float my_ALPHA,my_BETA;


    ChunkTask_GEN(int firstRow,int chunkSize,int rows,int cols,float *A,float *X,float *Y,float alpha,float beta,int incX,int incY) :
      my_firstRow(firstRow),my_CHUNKSIZE(chunkSize), my_ROWS(rows), my_COLS(cols), my_A(A), my_X(X), my_Y(Y), my_ALPHA(alpha), my_BETA(beta), my_INCX(incX), my_INCY(incY)
    {}
    task* execute()
//    void execute()
    {
        int cols = my_COLS;
        int rows = my_ROWS;
        float *a = my_A;

        int firstRow=my_firstRow;
        int lastRow=firstRow+my_CHUNKSIZE;
        int indexA=firstRow*cols;


        float *x = my_X;
        int incX=my_INCX;
        float alpha=my_ALPHA;

        int firstIndexX=0;
        if (incX < 0)
            firstIndexX=(cols-1)*incX*-1;
        int indexX=0;


        float *y = my_Y;
        int incY=my_INCY;
        float beta=my_BETA;

        int indexY=0;
        if (incY < 0)
            indexY=(rows-1)*incY*-1;
        indexY=indexY+firstRow*incY;

        int value=0;

        for ( int i=firstRow; i!=lastRow; ++i, indexY+=incY)
        {
#if REDUCE == 1
            ProductSum ps(a,x,indexA,firstIndexX,incX);
  #if PARTITIONER == 0
            parallel_reduce(blocked_range<size_t>(0,cols,GRAINSIZE), ps );
  #elif PARTITIONER == 1
            parallel_reduce(blocked_range<size_t>(0,cols,GRAINSIZE), ps,simple_partitioner());
  #else /* PARTITIONER is affinity*/
            static affinity_partitioner ap;
            parallel_reduce(blocked_range<size_t>(0,cols,GRAINSIZE), ps,ap );
  #endif
            y[indexY]=alpha*ps.my_value + beta*y[indexY];
            indexA+=cols;

#else
            value=0;
            indexX=firstIndexX;
            for (int j=0;j<cols;j++,indexX+=incX)
            {
                value+=a[indexA]*x[indexX];
                indexA=indexA+1;
            }
            y[indexY]=alpha*value + beta*y[indexY];
#endif
        }
        return NULL;
    }
};

class RootTask: public task
//class RootTask
{
public:
    int my_ROWS,my_COLS,my_KX,my_KY,my_INCX,my_INCY;
    int my_FACTOR;
    float my_ALPHA,my_BETA;
    float *my_A,*my_X,*my_Y;

    RootTask(int rows,int cols,float *A,float *X,float *Y,float alpha,float beta,int incX,int incY, int factor=1) :
                  my_ROWS(rows), my_COLS(cols), my_A(A), my_X(X), my_Y(Y), my_ALPHA(alpha), my_BETA(beta), my_INCX(incX), my_INCY(incY), my_FACTOR(factor)
    {}
    task* execute()
//    void execute()
    {
// Overrides virtual function task::execute
        long totalElements=my_ROWS * my_COLS;
        int cutoff=CUTOFF/my_FACTOR;
        int chunkCount=min((int)(totalElements/cutoff), PROCESSORS);
        if(chunkCount == 0)
          chunkCount=1;

#if SHOW_THREAD_INFO == 1
        cout<<endl<<"Dividing into "<<chunkCount<<" chunks"<<endl;
#endif
        ChunkTask_GEN* tasks[chunkCount];
        int a=my_COLS;
        int rows_remaining=my_ROWS;
        int i,chunkSize;
        for(i=0;i<chunkCount-1 ;i++)
        {
          long elements=totalElements/(chunkCount-i);

          if(elements < cutoff)
            break;

          chunkSize=elements/a;
          if(elements%a!=0)
            chunkSize=chunkSize+1;


          tasks[i] = new( allocate_child() ) ChunkTask_GEN(my_ROWS-rows_remaining,chunkSize,my_ROWS,my_COLS,my_A,my_X,my_Y,my_ALPHA,my_BETA,my_INCX,my_INCY);
//          tasks[i] = new ChunkTask_GEN(my_ROWS-rows_remaining,chunkSize,my_ROWS,my_COLS,my_A,my_X,my_Y,my_ALPHA,my_BETA,my_INCX,my_INCY);


          totalElements=totalElements - (chunkSize*a);

#if SHOW_THREAD_INFO == 1
        long operated_elements=(chunkSize*a);
        cout<<"\n Range: "<<my_ROWS-rows_remaining<<" to "<<my_ROWS-rows_remaining+chunkSize-1<<", with valid elements: "<<operated_elements<<", thead-id: "<<this_tbb_thread::get_id();
#endif

          rows_remaining-=chunkSize;
        }

        tasks[i] = new( allocate_child() ) ChunkTask_GEN(my_ROWS-rows_remaining,rows_remaining,my_ROWS,my_COLS,my_A,my_X,my_Y,my_ALPHA,my_BETA,my_INCX,my_INCY);
#if SHOW_THREAD_INFO == 1
        long operated_elements=(rows_remaining*a);
        cout<<"\n Range: "<<my_ROWS-rows_remaining<<" to "<<my_ROWS-rows_remaining+rows_remaining-1<<", with valid elements: "<<operated_elements<<", thead-id: "<<this_tbb_thread::get_id();
#endif
//        tasks[i] = new ChunkTask_GEN(my_ROWS-rows_remaining,rows_remaining,my_ROWS,my_COLS,my_A,my_X,my_Y,my_ALPHA,my_BETA,my_INCX,my_INCY);

        // Set ref_count to "total children plus one for the wait".
        set_ref_count(i+2);
        int j;

#if INNER_TIME == 1
        tick_count start,end,end1;
        start =tick_count::now();
#endif

          for(j=0;j<i;j++)
            spawn(*tasks[j]);
//            (tasks[j]->execute());    // Start tasks[j]

          // Start last task and wait for all children.
          spawn_and_wait_for_all(*tasks[j]);
#if INNER_TIME == 1
          end=tick_count::now();
          cout<<endl<<"Starting and finishing Chunk tasks:"<<(end-start).seconds()<<endl;
#endif
//          tasks[j]->execute();


        return NULL;
    }
};


class ProductSum
{
    float* my_a, *my_x;
    int my_baseA, my_baseX, my_incX, my_incA;
public:
    float my_value;
    void operator()( const blocked_range<size_t>& r )
    {
#if SHOW_INNER_THREAD_INFO == 1
        cout<<"\n range: "<<r.begin()<<" to "<<r.end()<<", object: "<<this<<", thead-id: "<<this_tbb_thread::get_id();
#endif
        float *a = my_a;
        float *x = my_x;
        int incX = my_incX;
        int incA = my_incA;
        int indexX = my_baseX;
        int indexA = my_baseA;

        float value = my_value;
        size_t end = r.end();

        for ( size_t i=r.begin(); i!=end; ++i,indexA+=incA,indexX+=incX )
            value += a[indexA]*x[indexX];

        my_value = value;
    }

    ProductSum( ProductSum& ps, split ) :
            my_a(ps.my_a), my_x(ps.my_x), my_baseA(ps.my_baseA), my_baseX(ps.my_baseX), my_incX(ps.my_incX), my_incA(ps.my_incA), my_value(0)
    {}
    void join( const ProductSum& ps )
    {
        my_value+=ps.my_value;
    }
    ProductSum(float *a, float *x,int baseA, int baseX, int incX, int incA=1 ) :
            my_a(a), my_x(x), my_baseA(baseA), my_baseX(baseX), my_incX(incX), my_incA(incA), my_value(0)
    {}
};

class TBB_SGEMV
{
    char my_TRANS;
    int my_ROWS,my_COLS,my_KX,my_KY,my_INCX,my_INCY;
    float my_ALPHA,my_BETA;
    float *my_A,*my_X,*my_Y;
public:
    void operator()( const blocked_range<size_t>& r ) const
    {
#if SHOW_THREAD_INFO == 1
        cout<<"\n range: "<<r.begin()<<" to "<<r.end()<<", object: "<<this<<", thead-id: "<<this_tbb_thread::get_id();
#endif
        int cols = my_COLS;
        int rows = my_ROWS;
        float *a = my_A;

        size_t indexA=r.begin()*cols;


        float *x = my_X;
        int incX=my_INCX;
        float alpha=my_ALPHA;

        int firstIndexX=0;
        if (incX < 0)
            firstIndexX=(cols-1)*incX*-1;
        int indexX=0;


        float *y = my_Y;
        int incY=my_INCY;
        float beta=my_BETA;

        int indexY=0;
        if (incY < 0)
            indexY=(rows-1)*incY*-1;
        indexY=indexY+r.begin()*incY;


        float value=0;




        for ( size_t i=r.begin(); i!=r.end(); ++i, indexY+=incY)
        {

#if REDUCE == 1
            ProductSum ps(a,x,indexA,firstIndexX,incX);
  #if PARTITIONER == 0
            parallel_reduce(blocked_range<size_t>(0,cols,GRAINSIZE), ps );
  #elif PARTITIONER == 1
            parallel_reduce(blocked_range<size_t>(0,cols,GRAINSIZE), ps,simple_partitioner());
  #else /* PARTITIONER is affinity*/
            static affinity_partitioner ap;
            parallel_reduce(blocked_range<size_t>(0,cols,GRAINSIZE), ps,ap );
  #endif
            y[indexY]=alpha*ps.my_value + beta*y[indexY];
            indexA+=cols;

#else
            indexX=firstIndexX;
            for (int j=0;j<cols;j++,indexX+=incX)
            {
                value+=a[indexA]*x[indexX];
                indexA=indexA+1;
            }
            y[indexY]=alpha*value + beta*y[indexY];
            value=0;
#endif
        }

    }

    TBB_SGEMV(int rows,int cols,float *A,float *X,float *Y,float alpha,float beta,int incX,int incY) :  my_ROWS(rows), my_COLS(cols), my_A(A), my_X(X), my_Y(Y), my_ALPHA(alpha), my_BETA(beta), my_INCX(incX), my_INCY(incY)
    {
//      printf("\n new SGEMV object: %p",(void *)this);
    }
};

class TBB_SGEMV_TRNS
{
    char my_TRANS;
    int my_ROWS,my_COLS,my_KX,my_KY,my_INCX,my_INCY;
    float my_ALPHA,my_BETA;
    float *my_A,*my_X,*my_Y;
public:
    void operator()( const blocked_range<size_t>& c ) const
    {
#if SHOW_THREAD_INFO == 1
        cout<<"\n TRNS range: "<<c.begin()<<" to "<<c.end()<<", object: "<<this<<", thead-id: "<<this_tbb_thread::get_id();
#endif
        int rows = my_ROWS;
        int cols = my_COLS;
        float *a = my_A;

        float *x = my_X;
        int incX=my_INCX;
        float alpha=my_ALPHA;

        int firstIndexX=0;
        if (incX < 0)
            firstIndexX=(cols-1)*incX*-1;
        firstIndexX=firstIndexX+c.begin()*incX;
        int indexX=firstIndexX;


        float *y = my_Y;
        int incY=my_INCY;
        float beta=my_BETA;

        int firstIndexY=0;
        if (incY < 0)
            firstIndexY=(rows-1)*incY*-1;
        int indexY=firstIndexY;

        size_t indexA=c.begin();
        size_t incrA=cols-c.end()+c.begin();

#if USE_TEMP_ARR == 0 /*Use temporary array*/
        int values_size=c.end()-c.begin();
        float *values=NULL;
        try
        {
            values=new float[values_size];
            int indexV=0;
            for (indexA=c.begin(), indexX=firstIndexX ; indexA!=c.end(); indexA++,indexV++)
            {
                values[indexV]=a[indexA]*x[indexX];
            }

            indexX+=incX;
            for (int i=1;i<rows;i++,indexX+=incX)
            {
                indexA+=incrA;

                for ( indexV=0;indexV<values_size;indexA++,indexV++)
                {
                    values[indexV]+=a[indexA]*x[indexX];
                }
            }

            indexY=firstIndexY;
            for ( indexV=0;indexV<values_size;indexV++, indexY+=incY)
            {
                y[indexY]=alpha*values[indexV] + beta*y[indexY];
            }
        }
        catch (...)
        {
            if (values!=NULL) delete []values;
            throw;
        }
        if (values!=NULL) delete []values;

#else
        float value=0;
        indexY=firstIndexY;
        for (size_t j=c.begin();j!=c.end();j++,indexY+=incY)
        {
  #if REDUCE == 1
            ProductSum ps(a,x,j,firstIndexX,incX,cols);
    #if PARTITIONER == 0
            parallel_reduce(blocked_range<size_t>(0,cols,GRAINSIZE), ps );
    #elif PARTITIONER == 1
            parallel_reduce(blocked_range<size_t>(0,cols,GRAINSIZE), ps,simple_partitioner() );
    #else /* PARTITIONER is affinity*/
            static affinity_partitioner ap;
            parallel_reduce(blocked_range<size_t>(0,cols,GRAINSIZE), ps,ap );
    #endif
            y[indexY]=alpha*ps.my_value + beta*y[indexY];

  #else
            indexX=firstIndexX;
            indexA=j;
            for ( int i=0; i<rows; ++i, indexX+=incX)
            {
                value+=a[indexA]*x[indexX];
                indexA=indexA+cols;
            }
            y[indexY]=alpha*value + beta*y[indexY];
            value=0;
  #endif
        }
#endif


    }

    TBB_SGEMV_TRNS(int rows,int cols,float *A,float *X,float *Y,float alpha,float beta,int incX,int incY) :  my_ROWS(rows), my_COLS(cols), my_A(A), my_X(X), my_Y(Y), my_ALPHA(alpha), my_BETA(beta), my_INCX(incX), my_INCY(incY)
    {
//      printf("\n new SGEMV_TRNS object: %p",(void *)this);
    }
};

void SGEMV(char TRANS,int M,int N,float ALPHA,float* A,int LDA,float* X,int INCX,float BETA,float* Y,int INCY)
{
#if USE_TASKS == 1
    if (TRANS=='N'||TRANS=='n')
    {
        RootTask& a = *new(task::allocate_root()) RootTask(M,N,A,X,Y,ALPHA,BETA,INCX,INCY,1);
//          RootTask& a = *new RootTask(N,N,A,X,Y,ALPHA,BETA,INCX,INCY);
        task::spawn_root_and_wait(a);
//          a.execute();
    }
#else

  #if PARTITIONER == 0
      if (TRANS=='N'||TRANS=='n')
      {
          parallel_for(blocked_range<size_t>(0,M,GRAINSIZE), TBB_SGEMV(M,N,A,X,Y,ALPHA,BETA,INCX,INCY));
      }
      else
      {
          parallel_for(blocked_range<size_t>(0,N,GRAINSIZE), TBB_SGEMV_TRNS(M,N,A,X,Y,ALPHA,BETA,INCX,INCY));
      }
  #elif PARTITIONER == 1
      if (TRANS=='N'||TRANS=='n')
      {
          parallel_for(blocked_range<size_t>(0,M,GRAINSIZE), TBB_SGEMV(M,N,A,X,Y,ALPHA,BETA,INCX,INCY),simple_partitioner());
      }
      else
      {
          parallel_for(blocked_range<size_t>(0,N,GRAINSIZE), TBB_SGEMV_TRNS(M,N,A,X,Y,ALPHA,BETA,INCX,INCY),simple_partitioner());
      }
  #else /* PARTITIONER is affinity*/
      static affinity_partitioner ap;

      if (TRANS=='N'||TRANS=='n')
      {
          parallel_for(blocked_range<size_t>(0,M,GRAINSIZE), TBB_SGEMV(M,N,A,X,Y,ALPHA,BETA,INCX,INCY),ap);
      }
      else
      {
          parallel_for(blocked_range<size_t>(0,N,GRAINSIZE), TBB_SGEMV_TRNS(M,N,A,X,Y,ALPHA,BETA,INCX,INCY),ap);
      }
  #endif

#endif



    return;
}



