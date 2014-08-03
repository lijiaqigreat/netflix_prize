#include "nsvd3.h"

#define MAX_RANK 5
#define ALPHA 0.5
#define BETA 0.25
#define DMAX 1e80
#define DMAXL 184
#define SWAP(a,b,tmp) tmp=a;a=b;b=tmp;
#define ABS(a) (a)>0?(a):-(a)
#define DIM 10

#define LAMDA 0.01

typedef double Double;
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

  }
  //update LAMDA
  int tt;
  add(dv,dv,v1,LAMDA,data.ndim);
  for(tt=0;tt<data.ndim;tt++){
    ddv[tt*data.ndim+tt]+=LAMDA;
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
  f=LAMDA*dot(v1,v1,ndim)/2;
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
int train_prepare_param(void** param_,int argc,char **argv)
{
  if(argc!=3){
    printf("train usage:\n dimension\n k1\n k2\n");
    exit(1);

  }
  Parameter* param=MALLOC(Parameter,1);
  *param_=param;
  int dim=atoi(argv[0]);
  Double k1,k2;
  sscanf(argv[1],"%lf",&k1);
  sscanf(argv[2],"%lf",&k2);
  param->dim=dim;
  param->k1=k1;
  param->k2=k2;
  param->line_alpha=0.5;
  param->line_beta=0.35;
  return sizeof(Parameter);
}

int train_prepare_vector(void **vector,void *param_)
{
  const Parameter* param=(const Parameter*) param_;
  int size=param->dim*2*(ncust+nmovie);
  Double *vector_=MALLOC(Double,size);
  *vector=vector_;
  int t;
  for(t=0;t<size;t++){
    vector_[t]=(rand()/(Double)RAND_MAX-0.5);
    if(t%(param->dim*2)>param->dim){
      vector_[t]/=1000;
    }

  }
  return size*sizeof(Double);
}

int train_prepare_report(void **report,void *param_)
{
  *report=MALLOC(Report_T,2);
  memset(*report,0,sizeof(Report_T)*2);
  return sizeof(Report_T)*2;
}

int train_prepare_data(void **datap,const Int64* data1,const void* param_)
{
  const Parameter* param=(const Parameter*)param_;
  void *datap_=malloc(sizeof(Int)*(ncust+nmovie+2)+sizeof(Int64)*(ntrain*2));
  *datap=datap_;
  int t1;
  memset(datap_,0,sizeof(Int)*(ncust+nmovie+2));
  Int* p1=(Int*)datap_;
  datap_+=sizeof(Int)*(ncust+1);
  Int* p2=(Int*)datap_;
  datap_+=sizeof(Int)*(nmovie+1);

  Int64* sort1=(Int64*)datap_;
  datap_+=sizeof(Int64)*(ntrain);
  Int64* sort2=(Int64*)datap_;
  datap_+=sizeof(Int64)*(ntrain);
  for(t1=0;t1<ntrain;t1++){
    Int64 n=data1[t1];
    Int n_movie=n&0xffff;
    Int n_cust=(n>>16)&0xffffff;
    //Int n_ranking=(n>>40)&0xff;
    //Int n_time=(n>>48)&0xffff;
    //collaberate with calculating accumulated sum
    p1[n_cust+1]++;
    p2[n_movie+1]++;
  }

  //calculate accumulated sum
  for(t1=1;t1<=ncust;t1++){
    p1[t1]+=p1[t1-1];
  }
  for(t1=1;t1<=nmovie;t1++){
    p2[t1]+=p2[t1-1];
  }
  //calculate sort
  for(t1=0;t1<ntrain;t1++){
    Int64 n=data1[t1];
    Int n_movie=n&0xffff;
    Int n_cust=(n>>16)&0xffffff;
    Int n_ranking=(n>>40)&0xff;
    Int n_time=(n>>48)&0xffff;
    sort1[p1[n_cust]++]=((Int64)n_movie<<32)+((Int64)n_ranking<<16)+n_time;
    sort2[p2[n_movie]++]=((Int64)n_cust<<32)+((Int64)n_ranking<<16)+n_time;
  }

  //recover p1,p2
  for(t1=ncust;t1>0;t1--){
    p1[t1]=p1[t1-1];
  }
  p1[0]=0;
  for(t1=nmovie;t1>0;t1--){
    p2[t1]=p2[t1-1];
  }
  p2[0]=0;
  if(p1[ncust]!=ntrain||p2[nmovie]!=ntrain){
    printf("error in preparing data: %d %d\n",p1[ncust],p2[nmovie]);
    exit(1);
  }
  return sizeof(Int)*(ncust+nmovie+2)+sizeof(Int64)*ntrain*2;
}

void train_print_param(char* rt,const void *param_)
{
  const Parameter* param=(const Parameter*)param_;
  rt+=sprintf(rt,"dimension: %d\n",param->dim);
  rt+=sprintf(rt,"k1,k2: %.4le,%.4le\n",param->k1,param->k2);
  rt+=sprintf(rt,"alpha,beta: %.4le,%.4le\n",param->line_alpha,param->line_beta);
}
void train_print_report(char* rt,const void *report)
{
  Report_T* reports=(Report_T*)report;
  int t;
  for(t=0;t<2;t++){
    rt+=sprintf(rt,"-- newton report %d --\n",t);
    rt+=sprintf(rt,"line_search,non_stable,singular: %d,%d,%d\n",reports[t].line_search,reports[t].non_stable,reports[t].singular);
    rt+=sprintf(rt,"size of dv,dvt: %le %le\n",reports[t].dv,reports[t].dvt);
    rt+=sprintf(rt,"dS, S1: %le %le\n",reports[t].dS,reports[t].S1);
    rt+=sprintf(rt,"time,training score: %lf,%lf\n",reports[t].time,reports[t].score);
  }
}

int train_prepare_vector(void** vector_,void* param)
{
  int size=(ncust+nmovie)*DIM;
  Double* vector=malloc(sizeof(Double)*size);
  *vector_=vector;
  for(size--;size>=0;size--){
    ((Double*)vector)[size]=0;
  }
  return size*sizeof(Double);
}

void train(const void* dat,const void* param,void *vector,void* report_x_2)
{
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

  

  //init vector
  Double* v1=(Double*)vector;
  vector+=sizeof(Double)*data.cust.n*data.ndim;
  Double* v2=(Double*)vector;
  vector+=sizeof(Double)*data.cust.n*data.ndim;
  if(v1[0]==0){
    //init vector
    int t1,t2;
    for(t1=0;t1<ncust;t1++){
      v1[t1*DIM]=1;
      for(t2=1;t2<DIM;t2++){
        v1[t1*DIM+t2]=((Double)rand()/RAND_MAX-0.5)/100;
      }
    }
    for(t1=0;t1<nmovie;t1++){
      Double sum=0;
      for(t2=data.movie.p[t1];t2<data.movie.p[t1+1];t2++){
        int rank=data.movie.i[t2]%5;
        sum+=rank;
      }
      sum/=data.movie.p[t1+1]-data.movie.p[t1];
      sum/=4;
      if(sum>0.99999){
        sum=DMAXL;
      }else if(sum<0.00001){
        sum=-DMAXL;
      }else{
        sum=-log((1/sum)-1)/2;
      }
      v2[t1*DIM]=sum;
      for(t2=1;t2<DIM;t2++){
        v2[t1*DIM+t2]=((Double)rand()/RAND_MAX-0.5)/100;
      }
    }
    return;
  }

  printf("starting iteration\n");
  time_t time0,time1;
  time(&time0);
  newton(v1,v2,data,data.cust,data.movie);
  newton(v2,v1,data,data.movie,data.cust);
  time(&time1);
  printf("time elapsed: %lu\n",time1-time0);
}


Double train_get_score(const void* vector,const void *param,int custi,int moviei,int time)
{
  Double* v1=(Double*)vector;
  vector+=sizeof(Double)*ncust*DIM;
  Double* v2=(Double*)vector;
  vector+=sizeof(Double)*nmovie*DIM;

  Double xy=dot(v1+custi*DIM,v2+moviei*DIM,DIM);
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
