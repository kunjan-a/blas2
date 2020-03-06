#include "blas2.h"
#include "tbb/tbb.h"
#include <cstdio>
#include <iostream>
#include <cmath>
using namespace std;
using namespace tbb;



class ChunkTask_SYM: public task
//class ChunkTask_SYM
{
public:
    int my_firstRow,my_ROWS,my_COLS,my_KX,my_KY,my_INCX,values_size;
    float *my_A,*my_X;
    float *values;


    ChunkTask_SYM(int firstRow,int rows,int cols,float *A,float *X,int incX) :
      my_firstRow(firstRow),my_ROWS(rows), my_COLS(cols), my_A(A), my_X(X), my_INCX(incX), values(NULL)
    {}
    task* execute()
//    void execute()
    {
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

        values_size=max(my_ROWS,cols-firstRow);
        try
        {
            int indexX;
            int indexA=firstRow*cols + firstRow;

            values=new float[values_size];
            for(int k=0;k<values_size;k++)
              values[k]=0;

            float value1_x;
            int indexV=0;
            int indexV1;
            float value=0;

            for (int i=firstRow;i!=lastRow;i++,indexV++)
            {
              indexX=firstIndexX;
              value1_x=x[indexX];
              value=a[indexA]*value1_x;
              indexA++;
              indexX+=incX;
              indexV1=indexV+1;

              for ( int j=i+1; j<cols; j++, indexA++,indexX+=incX,indexV1++)
              {
                  value+=a[indexA]*x[indexX];
                  values[indexV1]+=a[indexA]*value1_x;
              }

              values[indexV]+=value;
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

        return NULL;
    }
};

class SymmetricTask: public task
//class SymmetricTask
{
public:
    int my_ROWS,my_COLS,my_KX,my_KY,my_INCX,my_INCY;
    int my_FACTOR;
    float my_ALPHA,my_BETA;
    float *my_A,*my_X,*my_Y;

    SymmetricTask(int rows,int cols,float *A,float *X,float *Y,float alpha,float beta,int incX,int incY, int factor=1) :
                  my_ROWS(rows), my_COLS(cols), my_A(A), my_X(X), my_Y(Y), my_ALPHA(alpha), my_BETA(beta), my_INCX(incX), my_INCY(incY), my_FACTOR(factor)
    {}
    task* execute()
//    void execute()
    {
// Overrides virtual function task::execute
        long totalElements=((my_ROWS+1)*my_COLS)/2;
        int cutoff=CUTOFF/my_FACTOR;
        int chunkCount=min((int)(totalElements/cutoff), PROCESSORS);
        if(chunkCount == 0)
          chunkCount=1;

#if SHOW_THREAD_INFO == 1
        cout<<endl<<"Dividing into "<<chunkCount<<" chunks"<<endl;
#endif
        ChunkTask_SYM* tasks[chunkCount];
        int a=my_COLS;
        int rows_remaining=my_ROWS;
        int i,chunkSize;
        for(i=0;i<chunkCount-1 ;i++)
        {
          long elements=totalElements/(chunkCount-i);

          if(elements < cutoff)
            break;

          chunkSize=a - (double)(sqrt(a*(a+1) - 2*elements)-0.5);


          tasks[i] = new( allocate_child() ) ChunkTask_SYM(my_ROWS-rows_remaining,chunkSize,my_COLS,my_A,my_X,my_INCX);
//          tasks[i] = new ChunkTask_SYM(my_ROWS-rows_remaining,chunkSize,my_COLS,my_A,my_X,my_INCX);


          int l=a+1-chunkSize;
          totalElements=totalElements - (chunkSize*(a + l))/2;

#if SHOW_THREAD_INFO == 1
        long operated_elements=(chunkSize*(a + l))/2;
        cout<<"\n Range: "<<my_ROWS-rows_remaining<<" to "<<my_ROWS-rows_remaining+chunkSize-1<<", with valid elements: "<<operated_elements<<", thead-id: "<<this_tbb_thread::get_id();
#endif

          a=l-1;
          rows_remaining-=chunkSize;


        }

        tasks[i] = new( allocate_child() ) ChunkTask_SYM(my_ROWS-rows_remaining,rows_remaining,my_COLS,my_A,my_X,my_INCX);
#if SHOW_THREAD_INFO == 1
        long operated_elements=(rows_remaining*(1 + rows_remaining))/2;
        cout<<"\n Range: "<<my_ROWS-rows_remaining<<" to "<<my_ROWS-rows_remaining+rows_remaining-1<<", with valid elements: "<<operated_elements<<", thead-id: "<<this_tbb_thread::get_id();
#endif
//        tasks[i] = new ChunkTask_SYM(my_ROWS-rows_remaining,rows_remaining,my_COLS,my_A,my_X,my_INCX);

        // Set ref_count to "total children plus one for the wait".
        set_ref_count(i+2);
        int j;
        try
        {
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




          float *temp_values=new float[my_ROWS];
          int indexTemp=tasks[0]->my_firstRow;

          int values_size=tasks[0]->values_size;
          float *values=tasks[0]->values;
          int indexV;
          for(indexV=0;indexV<values_size;indexV++,indexTemp++){
            temp_values[indexTemp]=values[indexV];
          }



          for(j=1;j<=i;j++)
          {
            int values_size=tasks[j]->values_size;
            float *values=tasks[j]->values;
            indexTemp=tasks[j]->my_firstRow;

            for(int indexV=0;indexV<values_size;indexV++,indexTemp++){
              temp_values[indexTemp]+=values[indexV];
            }
          }


          int firstIndexY=0;
          if (my_INCY < 0)
              firstIndexY=(my_ROWS-1)*my_INCY*-1;
          int indexY=firstIndexY;

          for(indexV=0;indexV<my_ROWS;indexV++,indexY+=my_INCY)
          {
            my_Y[indexY]=my_BETA * my_Y[indexY] + my_ALPHA * temp_values[indexV];
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


void SSYMV(char UPLO,int N,float ALPHA,float* A,int LDA,float* X,int INCX,float BETA,float* Y,int INCY)
{
#if USE_TASKS == 1
    if (UPLO == 'u'||UPLO == 'U')
    {
        SymmetricTask& a = *new(task::allocate_root()) SymmetricTask(N,N,A,X,Y,ALPHA,BETA,INCX,INCY,2);
//          SymmetricTask& a = *new SymmetricTask(N,N,A,X,Y,ALPHA,BETA,INCX,INCY);
        task::spawn_root_and_wait(a);
//          a.execute();
    }
#else

        #if PARTITIONER == 0
//                parallel_for(blocked_range<size_t>(0,N,GRAINSIZE), TBB_SSYMV(N,N,A,X,Y,ALPHA,BETA,INCX,INCY));
        #elif PARTITIONER == 1
//                parallel_for(blocked_range<size_t>(0,N,GRAINSIZE), TBB_SSYMV(N,N,A,X,Y,ALPHA,BETA,INCX,INCY),simple_partitioner());
        #else
//            static affinity_partitioner ap;
//                parallel_for(blocked_range<size_t>(0,N,GRAINSIZE), TBB_SSYMV(N,N,A,X,Y,ALPHA,BETA,INCX,INCY),ap);
        #endif
#endif
        return;
    }



