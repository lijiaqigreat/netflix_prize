
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#define MAX_RANK 5
#define ALPHA 1.5
#define BETA 0.5 
#define DMAX 1e80
#define DMAXL 184
typedef double Double;
typedef uint32_t Int;
typedef struct{
  Int ndata;
  Int ncust;
  Int nmovie;
  Int ndim;
  Int *ind1;
  Int *ind2;
} Data;

//function g
static Double func_g_0(Double x);
//dot product
static Double dot(const Double *v1, const Double *v2,int ndim);
//get S given X, Y and constants
static Double getS_0(const Double *v1, const Double *v2, const Data* data);
static Double getScore(const Double *v1,const Double *v2, const Data* data);


static Double func_g_0(Double x)
{
  return x-log(cosh(x));
}

static Double dot(const Double *v1, const Double *v2,int ndim)
{
  Double f=0;
  int t1=0;
  for(t1=0;t1<ndim;t1++){
    f+=v1[t1]*v2[t1];
  }
  return f;
}

static void add(Double* v, const Double *v0, const Double *dv, Double step,int ndim){
  int t;
  for(t=0;t<ndim;t++){
    v[t]=v0[t]+dv[t]*step;
  }
}

static Double getSingleScore(const Double *v1,const Double *v2, const Data* data){
    Double xy=dot(v1,v2,data->ndim);
    Double score1;
    if(xy>DMAXL){
      score1=4;
    } else if(xy<-DMAXL){
      score1=0;
    } else{
      score1=4*exp(xy)/(exp(xy)+exp(-xy));
    }
    return score1;
}

static Double getScore(const Double *v1,const Double *v2, const Data* data)
{
  int t1;
  Double f=0;
  for(t1=0;t1<data->ndata;t1++){
    int custid=data->ind1[t1];
    int rankid=data->ind2[t1]%5;
    int movieid=(data->ind2[t1]/5);
    Double score1=getSingleScore(v1+custid*data->ndim,v2+movieid*data->ndim,data);
    Double error=score1-rankid;
    error*=error;
    f+=error;
  }
  return sqrt(f/data->ndata);
}

Double main_judge(void* qualify,void* vector){
    printf("starting judge\n");
    Data data;
    int t;

    //init constant
    do{
      Int* buffer=(Int*)qualify;
      qualify+=4*sizeof(Int);
      data.ndata=buffer[1];
      data.ncust=buffer[2];
      data.nmovie=buffer[3];
      data.ndim=10;

      data.ind1=qualify;
      qualify+=data.ndata*sizeof(Int);
      data.ind2=qualify;
    }while(0);



    /*
    //init accurence count
    Int* p1=(Int*)malloc(sizeof(Int)*nusers);
    Int* p2=(Int*)malloc(sizeof(Int)*nranks);
    memset(p1,0,sizeof(Int)*nusers);
    memset(p2,0,sizeof(Int)*nranks);
    int t;
    for(t=0;t<ndata;t++){
      p1[custind[t]]++;
      p2[rankind[t]]++;
    }
    */
    printf("ndata: %d\n",data.ndata);
    //data.ndata=10000000;


    //init vector
    Double* v1=(Double*)vector;
    Double* v2=v1+data.ncust*data.ndim;
    Double score=getScore(v1,v2,&data);
    printf("SCORE: %f\n",score);
}
