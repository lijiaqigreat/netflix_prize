#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

#define MAX_RANK 5
#define ALPHA 0.5
#define BETA 0.25
#define DMAX 1e80
#define DMAXL 184
#define MALLOC(type,size) (type*)malloc(sizeof(type)*(size));
#define SWAP(a,b,tmp) tmp=a;a=b;b=tmp;
#define ABS(a) (a)>0?(a):-(a)
#define DIM 10

#define LAMDA 0.00001

typedef double Double;
typedef uint32_t Int;
typedef struct{
  //number of cust/movie
  Int n;
  //accumulated accurence count
  Int *p;
  //corresponding movie/cust index and ranking
  Int *i;
} Index;
typedef struct{
  Int ndata;
  Index cust;
  Index movie;
  Int ndim;
} Data;

//function g,g',g''
static Double func_g_0(Double x);
static Double func_g_1(Double x);
static Double func_g_2(Double x);
//dot product
static Double dot(const Double *v1, const Double *v2,int ndim);
static void add(Double *v,const Double *v0,const Double *dv,Double step,int ndim);
//get S given X, Y and constants
static Double newton(Double *v1,const Double *v2,const Data data,const Index i1,const Index i2);
static void writeVec(const Double *v1,const Double *v2,Data data,const char* name);
static Double getS_0(const Double *v1, const Double *v2,const Int *i1,const Int *i2,int ndim);
static int updateV(Double S0,Double dots,Double *tmpv,Double *dv,Double *v1_,const Double *v2,const Int * i1,const Int *i2,int ndim);

void request_judge(void* vector);


static Double func_abs(Double x){
  return x<0?-x:x;
}

static Double func_g_0(Double x)
{
  if(x>DMAXL){
    return 0;
  }else if(x<-DMAXL){
    return 2*x;
  }else{
    return x-log(2*cosh(x));
  }
}

static Double func_g_1(Double x)
{
  return 1-tanh(x);
}

static Double func_g_2(Double x)
{
  Double tmp=cosh(x);
  return -1/(tmp*tmp);
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

static void add(Double* v, const Double *v0, const Double *dv, Double step,int ndim)
{
  int t;
  for(t=0;t<ndim;t++){
    v[t]=v0[t]+dv[t]*step;
  }
}

//print matrix
static void printm(const Double* m,int ndim)
{
  int t2;
  for(t2=0;t2<ndim*ndim;t2++){
    if(t2%ndim==0){
      printf("[");
    }
    printf("%e,",m[t2]);
    if(t2%ndim==(ndim-1)){
      printf("],\n");
    }
  }
}
static void printv(const Double* m,int ndim)
{
  int t2;
  for(t2=0;t2<ndim;t2++){
    printf("%e ",m[t2]);
  }
  printf("\n");
}

static int solve(Double *matrix,Double *b,const int ndim){
  int t1,t2,t3;
  int errors=0;
  for(t1=0;t1<ndim;t1++){
    //find max pivot 
    Double max=0;
    int maxi=-1;
    Double tmp;
    for(t2=t1;t2<ndim;t2++){
      tmp=func_abs(matrix[t2*ndim+t1]);
      if(tmp>max){
        maxi=t2;
        max=tmp;
      }
    }
    //swap to max
    for(t2=t1;t2<ndim;t2++){
      SWAP(matrix[t1*ndim+t2],matrix[maxi*ndim+t2],tmp);
    }
    SWAP(b[t1],b[maxi],tmp);
    if(max<1e-10){
      //fprintf(stderr,"LE warning!!!: %d %e\n",t1,max);
      //printm(matrix,ndim);
      //exit(1);
      errors+=1;
      matrix[t1*ndim+t1]=1;
      for(t2=t1+1;t2<ndim;t2++){
        matrix[t1*ndim+t2]=0;
      }
      b[t1]=0;
      continue;
    }
    //unify row
    tmp=matrix[t1*ndim+t1];
    for(t2=t1;t2<ndim;t2++){
      matrix[t1*ndim+t2]/=tmp;
    }
    b[t1]/=tmp;
    //reduce rows
    for(t2=0;t2<ndim;t2++){
      if(t2==t1){
        continue;
      }
      Double c=matrix[t2*ndim+t1];
      add(matrix+t2*ndim+t1,matrix+t2*ndim+t1,matrix+t1*ndim+t1,-c,ndim-t1);
      b[t2]-=c*b[t1];
    }
  }
  return errors;
  //printf("allgood\n");
}


/*
 * meaningful stats:
 * size
 * time
 * singluarity
 * deltav length
 * updateV info
 */
static int newton2(Double *fS,Double *v1,const Double *v2,const Data data,const Int *i1,const Int *i2){
  static int __tmp=-1;
  __tmp++;
  Double dv[DIM];
  Double deltav[DIM];
  Double tmpv[DIM];
  Double ddv[DIM*DIM];
  Double ddv2[DIM*DIM];
  memset(dv,0,sizeof(Double)*data.ndim);
  memset(ddv,0,sizeof(Double)*data.ndim*data.ndim);
  Double S0=0;
  Double fscore=0;
  const Int *i;
  //printf("size: %d\n",i1.p[t1+1]-i1.p[t1]);
  for(i=i1;i!=i2;i++){
    const Double *v2_=v2+(*i)/5*data.ndim;
    int rankid=(*i)%5;
    Double xy=dot(v1,v2_,data.ndim);
    //update score
    Double score=exp(func_g_0(xy));
    score=4*score-rankid;
    fscore+=score*score;
    //update S0
    S0-=rankid*func_g_0(xy)+(4-rankid)*func_g_0(-xy);
    Double c=-rankid*func_g_1(xy)+(4-rankid)*func_g_1(-xy);
    //update S1
    add(dv,dv,v2_,c,data.ndim);
    //update S2
    c=-rankid*func_g_2(xy)-(4-rankid)*func_g_2(-xy);
    int tt;
    for(tt=0;tt<data.ndim*data.ndim;tt++){
      ddv[tt]+=c*v2_[tt%data.ndim]*v2_[tt/data.ndim];
    }

    //update LAMDA
    add(dv,dv,v1,LAMDA,data.ndim);
    for(tt=0;tt<data.ndim;tt++){
      ddv[tt*data.ndim+tt]+=LAMDA;
    }
  }
  //solve equation
  //A*A
  int t2,t3;
  for(t2=0;t2<data.ndim*data.ndim;t2++){
    ddv2[t2]=0;
    for(t3=0;t3<data.ndim;t3++){
      ddv2[t2]+=ddv[t3*data.ndim+t2%data.ndim]*ddv[t2/data.ndim*data.ndim+t3];
    }
  }
  //backup dv
  memcpy(tmpv,dv,sizeof(Double)*data.ndim);
  //A*A*tmpv=dv
  solve(ddv2,tmpv,data.ndim);
  //deltav=-A*tmpv
  for(t2=0;t2<data.ndim;t2++){
    Double tmp=0;
    for(t3=0;t3<data.ndim;t3++){
      tmp-=ddv[t2*data.ndim+t3]*tmpv[t3];
    }
    deltav[t2]=tmp;
  }
  //dots<0
  Double dots=dot(deltav,dv,data.ndim);
  if(dots>0){
    fprintf(stderr,"error: %f %d %d %d\n",dots,(int)(i2-i1),__tmp,0);
    return 0;
  }
  *fS+=S0;
  return updateV(S0,dots,tmpv,deltav,v1,v2,i1,i2,data.ndim)-1;
}

static Double newton(Double *v1,const Double *v2,const Data data,const Index i1,const Index i2){
  Double fS=0;
  Double fscore=0;
  int tdata,t1,t2,t3;
  int vcount=0;

  for(t1=0;t1<i1.n;t1++){
    if(t1%10000==0){
      //printf("partial progress: %d %d\n",t1,vcount);
    }
    //get ddv,dv
    vcount+=newton2(&fS,v1+t1*data.ndim,v2,data,i1.i+i1.p[t1],i1.i+i1.p[t1+1]);
  }
  printf("S: %f\n",fS);
  printf("vcount: %d\n",vcount);
  return fS;
}
/*
 * meaningful stats:
 * iteration
 * delta S
 * new S
 */
static int updateV(Double S0,Double dots,Double *tmpv,Double *dv,Double *v1_,const Double *v2,const Int * i1,const Int *i2,int ndim){
  static int __tmp=-1;
  __tmp++;
  int count=0;
  Double step=1;
  memcpy(tmpv,v1_,sizeof(Double)*ndim);
  Double S1;
  while(1){
    count++;
    if(step<1e-100){
      int t;
      Double max=0;
      for(t=0;t<ndim;t++){
        if(func_abs(dv[t])>max){
          max=func_abs(dv[t]);
        }
      }
      fprintf(stderr,"error max: %d %e %e %e\n",__tmp,max,dots,S0-S1);
      memcpy(v1_,tmpv,sizeof(Double)*ndim);
      return count;
    }
    add(v1_,tmpv,dv,step,ndim);
    S1=getS_0(v1_,v2,i1,i2,ndim);
    if(S0==S1||S0-S1>=-step*dots*BETA){
      return count;
    }
    step*=ALPHA;
  }
}
static Double getS_0(const Double *v1, const Double *v2,const Int *i1,const Int *i2,int ndim)
{
  int t1;
  Double f;
  f=0;
  while(i1!=i2){
    int movieid=(*i1)/5;
    int rankid=(*i1)%5;
    Double xy=dot(v1,v2+movieid*ndim,ndim);
    Double c=-rankid*func_g_0(xy)-(4-rankid)*func_g_0(-xy);
    f+=c;
    i1++;
  }
  return f;
}
static void writeVec(const Double *v1,const Double *v2,Data data,const char* name){
  printf("write vector\n");
  FILE *file=fopen(name,"wb");
  fwrite(v1,sizeof(Double),data.cust.n*data.ndim,file);
  fwrite(v2,sizeof(Double),data.movie.n*data.ndim,file);
  fclose(file);
  printf("write vector done\n");
}

void main_train(int iter,void* dat,void *vector){
  printf("starting train\n");
  Data data;
  int t;

  //init constant
  do{
    Int* buffer=(Int*)dat;
    dat+=sizeof(Int)*3;
    data.ndata=buffer[0];
    data.cust.n=buffer[1];
    data.movie.n=buffer[2];
    data.ndim=10;
    data.cust.p=(Int*)dat;
    dat+=sizeof(Int)*(data.cust.n+1);
    data.movie.p=(Int*)dat;
    dat+=sizeof(Int)*(data.movie.n+1);
    data.cust.i=(Int*)dat;
    dat+=sizeof(Int)*data.ndata;
    data.movie.i=(Int*)dat;
    dat+=sizeof(Int)*data.ndata;
  }while(0);
  printf("(ndata: %d)(ncust: %d)(nmovie: %d)\n",data.ndata,data.cust.n,data.movie.n);
  //data.ndata=10000000;
  for(t=0;t<10;t++){
    printf("p2[%d]=%d\n",t,data.movie.p[t]);
  }

  if(vector==0){
    int size=(data.cust.n+data.movie.n)*data.ndim;
    vector=malloc(sizeof(Double)*size);
    for(size--;size>=0;size--){
      ((Double*)vector)[size]=((Double)rand()/RAND_MAX-0.5)/10;
    }
  }
  

  //init vector
  Double* v1=(Double*)vector;
  vector+=sizeof(Double)*data.cust.n*data.ndim;
  Double* v2=(Double*)vector;
  vector+=sizeof(Double)*data.cust.n*data.ndim;



  printf("starting iteration\n");
  time_t time0,time1;
  time(&time0);
  int iterleft=iter;
  while(iterleft!=0){
    newton(v1,v2,data,data.cust,data.movie);
    request_judge(v1);
    newton(v2,v1,data,data.movie,data.cust);
    request_judge(v1);
    printf("iterleft: %d\n",iterleft);
    time(&time1);
    printf("time elapsed: %lu\n",time1-time0);
    iterleft--;
  }
}
