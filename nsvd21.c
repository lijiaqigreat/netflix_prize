#include "nsvd2.h"
#include <sys/time.h>

#define SWAP(a,b,tmp) tmp=a;a=b;b=tmp;

typedef struct{
  int dim;
  Double k1;
  Double k2;
  Double line_alpha;
  Double line_beta;
} Parameter;

typedef struct{
  //updateV
  int line_search;
  //newton
  int non_stable;
  //newton
  int singular;
  Double dv;
  Double dvt;
  //updateV
  Double dS;
  //updateV
  Double S1;
  //newton
  Double score;
  //newton
  Double time;
} Report_T;

/*
 * TODO:
 * check dim for dim2
 * check each report element
 */
static inline Double func_abs(Double x);
static Double func_g_0(Double x);
static Double func_g_1(Double x);
static Double func_g_2(Double x);

static void add(Double *v,const Double *v0,const Double *dv,Double step,int ndim);
static Double dot(const Double *v1, const Double *v2,int ndim);
static Double dott(const Double *v1,const Double *v2,Int time,int dim);

static int solve(Double *matrix,Double *b,const int ndim);
static Double updateV(Report_T *report,Double *v1_,Double *tmpv,const Double* v2sort,const Double *dv,const Int64* i,const int N,const Double S0,const Double dots,const Parameter *param);
static Double getS_0(const Double* v1_,const Double* v2sort,const Int64 *i,const int N,int ndim);
static void newton(Report_T *report,Double* v1,const Double* v2,const int N,const Int* p,const Int64* i,const Parameter* param);

static inline Double func_abs(Double x)
{
  return x>0?x:-x;
}

static Double func_g_0(Double x)
{
  //log(double_max)
  if(x<-200.){
    return 2*x;
  }else{
    return -log(1+exp(-2*x));
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

static void add(Double* v, const Double *v0, const Double *dv, Double step,int ndim)
{
  int t;
  for(t=0;t<ndim;t++){
    v[t]=v0[t]+dv[t]*step;
  }
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

static Double dott(const Double *v1,const Double *v2,Int time,int dim)
{
  Double f=0;
  int t;
  for(t=0;t<dim;t++){
    f+= (v1[t]+v1[t+dim]*t)*
        (v2[t]+v2[t+dim]*t);
  }
  return f;
}

static int solve(Double *matrix,Double *b,const int ndim)
{
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

      for(t3=t1;t3<ndim;t3++){
        matrix[t2*ndim+t3]-=c*matrix[t1*ndim+t3];
      }
      b[t2]-=c*b[t1];
    }
  }
  return errors;
}

static Double updateV(Report_T *report,Double *v1_,Double *tmpv,const Double* v2sort,const Double *dv,const Int64* i,const int N,const Double S0,const Double dots,const Parameter *param)
{
  int dim2=param->dim*2;
  int count=0;
  Double step=1;
  memcpy(tmpv,v1_,sizeof(Double)*dim2);
  Double S1;
  while(1){
    add(tmpv,v1_,dv,step,dim2);
    S1=getS_0(tmpv,v2sort,i,N,dim2/2);
    printf("S0,dS,ldS,pdedS: %f %e %e %e\n",S0,S0-S1,step*dots,(S0-S1+dots*step*(2-step)/2)/((S0-S1)*(S0-S1)));
    if(S0-S1>=-step*dots*param->line_beta){
      memcpy(v1_,tmpv,sizeof(Double)*dim2);
      report->line_search+=count;
      report->dS+=S1-S0;
      report->S1+=S1;
      return step;
    }
    count++;
    step*=param->line_alpha;
  }
}

static Double getS_0(const Double* v1_,const Double* v2sort,const Int64 *i,const int N,int dim)
{
  int t;
  Double f=0;
  for(t=0;t<N;t++){
    Int64 n=i[t];
    //not needed since v2 is sorted
    //Int n_i=(n>>32)&0xffffff;
    Int n_ranking=(n>>16)&0xffff;
    Int n_time=n&0xffff;
    Double xy=dott(v1_,v2sort+t*dim*2,n_time,dim);
    //f-=rankid*func_g_0(xy)+(4-rankid)*func_g_0(-xy);
    f-=4*func_g_0(xy)-(8-2*n_ranking)*xy;
  }
  return f;
}

static void newton(Report_T *report,Double* v1,const Double* v2,const int N,const Int* p,const Int64* i,const Parameter* param)
{

  struct timeval time_start,time_end;
  gettimeofday(&time_start,0);
  int dim=param->dim;
  int dim2=dim*2;
  Double *dv=MALLOC(Double,dim2);
  Double *deltav=MALLOC(Double,dim2);
  Double *tmpv=MALLOC(Double,dim2);
  Double *ddv=MALLOC(Double,dim2*dim2);
  Double *ddv2=MALLOC(Double,dim2*dim2);
  int t1,t2,t3;
  //loop v1
  for(t1=0;t1<N;t1++){
    do{
      char tmpc[1000];
      //TODO remove
      //train_print_report(tmpc,report);
      //printf("t1: %d\n%s",t1,tmpc);
      if(t1==10){
        exit(0);
      }
    }while(0);
    Double S0=0;
    Double fscore=0;
    Double *v1_=v1+dim2*t1;
    memset(dv,0,sizeof(Double)*dim2);
    memset(ddv,0,sizeof(Double)*dim2*dim2);
    Double *v2sort=MALLOC(Double,dim2*(p[t1+1]-p[t1]));
    //init v2sort
    //TODO condense v2sort from dim2 to dim
    for(t2=p[t1];t2<p[t1+1];t2++){
      Int64 n=i[t2];
      Int n_i=(n>>32)&0xffffff;
      //Int n_ranking=(n>>16)&0xffff;
      //Int n_time=n&0xffff;
      memcpy(v2sort,v2+n_i*dim2,sizeof(Double)*dim2);
      v2sort+=dim2;
    }
    //restore v2sort
    v2sort-=dim2*(p[t1+1]-p[t1]);

    for(t2=0;t2<p[t1+1]-p[t1];t2++){
      Int64 n=i[t2+p[t1]];
      //Int n_i=(n>>32)&0xffffff;
      Int n_ranking=(n>>16)&0xffff;
      Int n_time=n&0xffff;
      const Double *v2_=v2sort+t2*dim2;
      
      Double xy=dott(v1_,v2_,n_time,dim);
      //update score
      Double score=4/(1+exp(-xy*2))-n_ranking;
      report->score+=score*score;
      //simplify all results
      //update S0
      //S0-=rankid*func_g_0(xy)+(4-rankid)*func_g_0(-xy);
      S0-=4*func_g_0(xy)-(8-2*n_ranking)*xy;
      //Double c=-rankid*func_g_1(xy)+(4-rankid)*func_g_1(-xy);
      Double c=4*tanh(xy)+((int)4-2*n_ranking);
      //update S1
      for(t3=0;t3<dim;t3++){
        tmpv[t3]=v2_[t3]+v2_[t3+dim]*n_time;
        dv[t3]+=c*tmpv[t3];
        dv[t3+dim]+=c*tmpv[t3]*n_time;
      }
      //update S2
      //c=-rankid*func_g_2(xy)-(4-rankid)*func_g_2(-xy);
      c=-4*func_g_2(xy);

      for(t3=0;t3<dim*dim;t3++){
        Double tmp=c*tmpv[t3%dim]*tmpv[t3/dim];
        ddv[t3]+=tmp;
        ddv[t3+dim]+=tmp*n_time;
        ddv[t3+dim*dim2]+=tmp*n_time;
        ddv[t3+dim*(dim2+1)]+=tmp*n_time*n_time;
      }
    }
    //CHANGE move this one loop out
    //update 
    Double k1=param->k1;
    Double k2=param->k1;
    for(t2=0;t2<dim;t2++){
      dv[t2]+=v1_[t2]*k1;
      ddv[t2*dim2+t2]+=k1;
      dv[t2+dim]+=v1_[t2+dim]*k2;
      ddv[(t2+dim)*(dim2+1)]+=k2;
    }
    //solve equation
    //A*A
    for(t2=0;t2<dim2*dim2;t2++){
      ddv2[t2]=0;
      for(t3=0;t3<dim2;t3++){
        ddv2[t2]+=ddv[t3*dim2+t2%dim2]*ddv[t2/dim2*dim2+t3];
      }
    }
    //backup dv
    memcpy(tmpv,dv,sizeof(Double)*dim2);
    //A*A*tmpv=dv
    report->singular+=solve(ddv2,tmpv,dim2);
    //deltav=-A*tmpv
    for(t2=0;t2<dim2;t2++){
      Double tmp=0;
      for(t3=0;t3<dim2;t3++){
        tmp-=ddv[t2*dim2+t3]*tmpv[t3];
      }
      deltav[t2]=tmp;
    }
    //dots<0
    Double dots=0;
    for(t2=0;t2<dim2;t2++){
      dots+=deltav[t2]*dv[t2];
    }
    if(dots>0){
      report->non_stable+=1;
      report->S1+=S0;
      free(v2sort);
      continue;
    }
    Double step = updateV(report,v1_,tmpv,v2sort,deltav,i+p[t1],p[t1+1]-p[t1],S0,dots,param);
    report->dv+=step*step*dot(deltav,deltav,dim);
    report->dvt+=step*step*dot(deltav+dim,deltav+dim,dim);
    free(v2sort);
    if(t1%1000==0){
      printf("debug score: %d %le %le\n",t1,report->S1,report->dS);
    }
  }
  free(dv);
  free(deltav);
  free(tmpv);
  free(ddv);
  free(ddv2);
  gettimeofday(&time_end,0);
  report->time=(time_end.tv_sec-time_start.tv_sec)+(time_end.tv_usec-time_end.tv_usec)*1e-6;
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

void train(const void* datap,const void* param_, void* vector_,void* report_x_2_)
{
  //convert pointer
  const Parameter* param=(const Parameter*)param_;
  Double* vector=(Double*)vector_;
  Report_T* report_x_2=(Report_T*)report_x_2_;

  //expand datap
  const void* datap_=datap;
  const Int* p1=(const Int*)datap_;
  datap_+=sizeof(Int)*(ncust+1);
  const Int* p2=(const Int*)datap_;
  datap_+=sizeof(Int)*(nmovie+1);

  const Int64* sort1=(const Int64*)datap_;
  datap_+=sizeof(Int64)*(ntrain);
  const Int64* sort2=(const Int64*)datap_;
  datap_+=sizeof(Int64)*(ntrain);

  Double* v1=vector;
  Double* v2=vector+param->dim*2*ncust;
  newton(report_x_2  ,v1,v2,ncust ,p1,sort1,param);
  newton(report_x_2+1,v2,v1,nmovie,p2,sort2,param);
}

Double train_get_score(const void* vector_,const void *param_,int custi,int moviei,int time)
{
  const Parameter* param=(const Parameter*)param_;
  int dim2=param->dim*2;
  Double* vector=(Double*)vector_;
  Double* v1=vector;
  Double* v2=vector+ncust*dim2;
  Double xy=dott(v1+custi*dim2,v2+moviei*dim2,time,dim2/2);
  return 4/(1+exp(-xy*2));
}
