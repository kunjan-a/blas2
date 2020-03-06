#include "blas2.h"
#include "tbb/tbb.h"
#include "tbb/tick_count.h"
#include <cstdio>
#include <iostream>
#include <cmath>
using namespace std;
using namespace tbb;

class ChunkTask_TRNG: public task
//class ChunkTask_TRNG
{
public:
    int my_firstRow,my_ROWS,my_COLS,my_KX,my_INCX,values_size;
    float *my_A,*my_X;
    float *values;
    tick_count my_start;


    ChunkTask_TRNG(int firstRow,int rows,int cols,float *A,float *X,int incX, tick_count start) :
      my_firstRow(firstRow),my_ROWS(rows), my_COLS(cols), my_A(A), my_X(X), my_INCX(incX), values(NULL),my_start(start)
    {}
    task* execute()
//    void execute()
    {
#if INNER_TIME == 1
        tick_count t1,t2;
        t1=tick_count::now();
#endif

        int firstRow=my_firstRow;
        int lastRow=firstRow+my_ROWS;

        int rows = my_COLS;
        int cols = my_COLS;
        float *a = my_A;

        float *x = my_X;
        int incX=my_INCX;

        int firstIndexX=0;
        if (incX < 0)
            firstIndexX=(cols-1)*incX*-1;
        firstIndexX=firstIndexX+firstRow*incX;

        if(firstRow==0)
        {
            int indexX;
            int indexA=firstRow*cols + firstRow;

//            values=new float[values_size];

            float value=0;

            for (int i=firstRow;i!=lastRow;i++)
            {
              indexX=firstIndexX;
              value=a[indexA]*x[indexX];
              indexA++;
              indexX+=incX;

              for ( int j=i+1; j<cols; j++, indexA++,indexX+=incX)
              {
                  value+=a[indexA]*x[indexX];
              }
              x[firstIndexX]=value;
              indexA+=i+1;
              firstIndexX+=incX;
            }
        }else
        {
          values_size=my_ROWS;
          try
          {
              int indexX;
              int indexA=firstRow*cols + firstRow;

              values=new float[values_size];

              int indexV=0;
              float value=0;
              for (int i=firstRow;i!=lastRow;i++,indexV++)
              {
                indexX=firstIndexX;
                value=a[indexA]*x[indexX];
                indexA++;
                indexX+=incX;

                for ( int j=i+1; j<cols; j++, indexA++,indexX+=incX)
                {
                    value+=a[indexA]*x[indexX];
                }

                values[indexV]=value;
                indexA+=i+1;
                firstIndexX+=incX;
              }
            }
            catch (...)
            {
                if (values!=NULL) delete []values;
                values=NULL;
                throw;
            }
        }
#if INNER_TIME == 1
        t2=tick_count::now();
        cout<<endl<<"Starting and finishing Chunk task:"<<(t1-my_start).seconds()<<" to "<<(t2 - my_start).seconds()<<" firstRow:"<<firstRow<<endl;
#endif

        return NULL;
    }
};

class TriangleTask: public task
//class TriangleTask
{
public:
    int my_ROWS,my_COLS,my_KX,my_INCX;
    int my_FACTOR;
    float *my_A,*my_X;

    TriangleTask(int rows,int cols,float *A,float *X,int incX, int factor=1) :
                  my_ROWS(rows), my_COLS(cols), my_A(A), my_X(X), my_INCX(incX), my_FACTOR(factor)
    {}
    task* execute()
//    void execute()
    {
// Overrides virtual function task::execute
        tick_count start;
#if INNER_TIME == 1
        tick_count start1,end,end1;
        start =tick_count::now();
#endif
        long totalElements=((my_ROWS+1)*my_COLS)/2;
        int cutoff=CUTOFF/my_FACTOR;
        int chunkCount=min((int)(totalElements/cutoff), PROCESSORS);
        if(chunkCount == 0)
          chunkCount=1;
#if SHOW_THREAD_INFO == 1
        cout<<endl<<"Dividing into "<<chunkCount<<" chunks"<<endl;
#endif
        ChunkTask_TRNG* tasks[chunkCount];
        int a=my_COLS;
        int rows_remaining=my_ROWS;
        int i,chunkSize;
        for(i=0;i<chunkCount-1 ;i++)
        {
          long elements=totalElements/(chunkCount-i);

          if(elements < cutoff)
            break;

          chunkSize=a - (double)(sqrt(a*(a+1) - 2*elements)-0.5);


          tasks[i] = new( allocate_child() ) ChunkTask_TRNG(my_ROWS-rows_remaining,chunkSize,my_COLS,my_A,my_X,my_INCX,start);
//          tasks[i] = new ChunkTask_TRNG(my_ROWS-rows_remaining,chunkSize,my_COLS,my_A,my_X,my_INCX);


          int l=a+1-chunkSize;
          totalElements=totalElements - (chunkSize*(a + l))/2;

#if SHOW_THREAD_INFO == 1
        long operated_elements=(chunkSize*(a + l))/2;
        cout<<"\n Range: "<<my_ROWS-rows_remaining<<" to "<<my_ROWS-rows_remaining+chunkSize-1<<", with valid elements: "<<operated_elements<<", thead-id: "<<this_tbb_thread::get_id();
#endif

          a=l-1;
          rows_remaining-=chunkSize;


        }

        tasks[i] = new( allocate_child() ) ChunkTask_TRNG(my_ROWS-rows_remaining,rows_remaining,my_COLS,my_A,my_X,my_INCX,start);
#if SHOW_THREAD_INFO == 1
        long operated_elements=(rows_remaining*(1 + rows_remaining))/2;
        cout<<"\n Range: "<<my_ROWS-rows_remaining<<" to "<<my_ROWS-rows_remaining+rows_remaining-1<<", with valid elements: "<<operated_elements<<", thead-id: "<<this_tbb_thread::get_id();
#endif
//        tasks[i] = new ChunkTask_TRNG(my_ROWS-rows_remaining,rows_remaining,my_COLS,my_A,my_X,my_INCX);

        // Set ref_count to "total children plus one for the wait".
        set_ref_count(i+2);
        int j;
        try
        {
#if INNER_TIME == 1
        start1 =tick_count::now();
#endif
          for(j=0;j<i;j++)
            spawn(*tasks[j]);
//            (tasks[j]->execute());    // Start tasks[j]

          // Start last task and wait for all children.
          spawn_and_wait_for_all(*tasks[j]);
#if INNER_TIME == 1
          end=tick_count::now();
          cout<<endl<<"Starting and finishing Chunk tasks:"<<(end-start1).seconds()<<endl;
#endif

//          tasks[j]->execute();




//          int values_size=tasks[0]->values_size;
//          float *values=tasks[0]->values;
//          int indexV;
//          int firstIndexX=0;
//          if (my_INCX < 0)
//              firstIndexX=(my_COLS-1)*my_INCX*-1;
//          int indexX=firstIndexX+tasks[0]->my_firstRow*my_INCX;
//
//          for(indexV=0;indexV<values_size;indexV++,indexX+=my_INCX){
//            my_X[indexX]=values[indexV];
//          }

          int values_size=0;
          float *values=NULL;
          int indexV=0;
          int firstIndexX=0;
          if (my_INCX < 0)
              firstIndexX=(my_COLS-1)*my_INCX*-1;
          int indexX=i<1?0:firstIndexX+tasks[1]->my_firstRow*my_INCX;

          for(j=1;j<=i;j++)
          {
            values_size=tasks[j]->values_size;
            values=tasks[j]->values;
//            if(indexX!=firstIndexX+tasks[j]->my_firstRow*my_INCX)
//              cout<<endl<<"Error for j="<<j<<". indexX="<<indexX<<" new value = "<<(firstIndexX+tasks[j]->my_firstRow*my_INCX)<<endl;
            indexX=firstIndexX+tasks[j]->my_firstRow*my_INCX;

            for(indexV=0;indexV<values_size;indexV++,indexX+=my_INCX){
              my_X[indexX]=values[indexV];
            }
          }
#if INNER_TIME == 1
          end1=tick_count::now();
          cout<<"Reduce operation:"<<(end1-end).seconds()<<endl;
#endif

        }catch(...)
        {
          for(j=0;j<=i;j++)
          {
            if(tasks[j]->values!=NULL)
            {
              delete []tasks[j]->values;
              tasks[j]->values=NULL;
            }
          }
        }

          for(j=0;j<=i;j++)
          {
            if(tasks[j]->values!=NULL)
            {
              delete []tasks[j]->values;
              tasks[j]->values=NULL;
            }
          }



        return NULL;
    }
};

void STRMV(char UPLO,char TRANS,char DIAG,int N,float*A,int LDA,float* X,int INCX)
{
#if USE_TASKS == 1
    if ((UPLO == 'u'||UPLO == 'U') && (TRANS == 'n' || TRANS == 'N') && (DIAG == 'n' || DIAG == 'N'))
    {
        TriangleTask& a = *new(task::allocate_root()) TriangleTask(N,N,A,X,INCX,1);
//          TriangleTask& a = *new TriangleTask(N,N,A,X,Y,ALPHA,BETA,INCX,INCY);
        task::spawn_root_and_wait(a);
//          a.execute();
    }
#else

        #if PARTITIONER == 0
//                parallel_for(blocked_range<size_t>(0,N,GRAINSIZE), TBB_STRMV(N,N,A,X,Y,ALPHA,BETA,INCX,INCY));
        #elif PARTITIONER == 1
//                parallel_for(blocked_range<size_t>(0,N,GRAINSIZE), TBB_STRMV(N,N,A,X,Y,ALPHA,BETA,INCX,INCY),simple_partitioner());
        #else
//            static affinity_partitioner ap;
//                parallel_for(blocked_range<size_t>(0,N,GRAINSIZE), TBB_STRMV(N,N,A,X,Y,ALPHA,BETA,INCX,INCY),ap);
        #endif
#endif

return;


}




