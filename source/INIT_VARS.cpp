#include "blas2.h"
#include "tbb/tbb.h"

using namespace tbb;

class TBB_INIT_A_GEN
{
    int my_COLS,my_ROWS;
    float my_MIN;
    float *my_A;
public:
    void operator()( const blocked_range<size_t>& r ) const
    {
  //        cout<<"\n range: "<<r.begin()<<" to "<<r.end()<<", object: "<<this<<", thead-id: "<<this_tbb_thread::get_id();
      int cols=my_COLS;
      int rows=my_ROWS;
      float minimum=my_MIN;
      float *A=my_A;

      size_t index=r.begin()*cols;

      for (int i=r.begin();i!=r.end();++i)
      {
          for (int j=0;j<cols;j++)
          {
              A[index]=rand()%VALUE_RANGE + minimum;
              index=index+1;
          }
      }
    }


    TBB_INIT_A_GEN(float *A, int rows, int cols, float minimum) :  my_A(A), my_ROWS(rows), my_COLS(cols), my_MIN(minimum)
    {
//      printf("\n new INIT_A_GEN object: %p",(void *)this);
    }
};

class TBB_INIT_A_SYM
{
    int my_DIM;
    float my_MIN;
    float *my_A;
public:
    void operator()( const blocked_range<size_t>& r ) const
    {
  //        cout<<"\n range: "<<r.begin()<<" to "<<r.end()<<", object: "<<this<<", thead-id: "<<this_tbb_thread::get_id();
      int cols=my_DIM;
      int rows=my_DIM;
      float minimum=my_MIN;
      float *A=my_A;

      size_t index=r.begin()*cols;
      size_t index2;

      for (int i=r.begin();i!=r.end();++i)
      {
          index+=i;
          index2=index+cols;
          A[index]=rand()%VALUE_RANGE + minimum;
          index++;
          for (int j=i+1;j<cols;j++)
          {
              A[index]=rand()%VALUE_RANGE + minimum;
              A[index2]=A[index];
              index=index+1;
              index2=index2+cols;
          }
      }
    }

    TBB_INIT_A_SYM(float *A, int dimension, float minimum) :  my_A(A), my_DIM(dimension), my_MIN(minimum)
    {
//      printf("\n new INIT_A_SYM object: %p",(void *)this);
    }
};

class TBB_INIT_A_TRNG
{
    int my_DIM;
    float my_MIN;
    float *my_A;
public:
    void operator()( const blocked_range<size_t>& r ) const
    {
  //        cout<<"\n range: "<<r.begin()<<" to "<<r.end()<<", object: "<<this<<", thead-id: "<<this_tbb_thread::get_id();
      int cols=my_DIM;
      float minimum=my_MIN;
      float *A=my_A;

      size_t index=r.begin()*cols;
      size_t index2;

      for (int i=r.begin();i!=r.end();++i)
      {
          index+=i;
          index2=index+cols;
          A[index]=rand()%VALUE_RANGE + minimum;
          for (int j=i+1;j<cols;j++)
          {
              A[++index]=rand()%VALUE_RANGE + minimum;
              A[index2]=0;
              index2+=cols;
          }
      }
    }

    TBB_INIT_A_TRNG(float *A, int dimension, float minimum) :  my_A(A), my_DIM(dimension), my_MIN(minimum)
    {
//      printf("\n new INIT_A_TRNG object: %p",(void *)this);
    }
};

void init_gen(float *A, int rows, int cols, float minimum)
{
  parallel_for(blocked_range<size_t>(0,rows), TBB_INIT_A_GEN(A,rows,cols,minimum));
}

void init_sym(float *A, int rows, float minimum)
{
  parallel_for(blocked_range<size_t>(0,rows), TBB_INIT_A_SYM(A,rows,minimum));
}

void init_trng(float *A, int rows, float minimum)
{
  parallel_for(blocked_range<size_t>(0,rows), TBB_INIT_A_TRNG(A,rows,minimum));
}
