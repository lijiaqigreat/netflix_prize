#include "nsvd2.h"

static Double main_judge(const Int64* data2, const void* vector, const void* param);

//data1
//data2
//dim
//iter
//k1
//k2
static void read_data(Int64** data,const char* path,int size);

static void read_data(Int64** data,const char* path,int size)
{
  FILE *file=fopen(path,"rb");
  *data=MALLOC(Int64,size);
  fread(*data,sizeof(Int64),size,file);
}

int main(int argc,char** argv)
{
  if(argc!=8){
    printf("usage:\n data1 path\n data2 path\n vector path\n iteration\n dimension\n k1\n k2\n");
  }
  
  char* tmpc=MALLOC(char,0x400);
  char* path1=argv[1];
  char* path2=argv[2];
  char* pathv=argv[3];
  int niter=atoi(argv[4]);
  Int64* data1;
  Int64* data2;
  read_data(&data1,path1,ntrain);
  read_data(&data2,path2,nprobe);

  void* param;
  void* vector;
  void* report;
  void* data;
  int param_s=train_prepare_param(&param,argc-5,argv+5);
  int vector_s=train_prepare_vector(&vector,param);
  int report_s=train_prepare_report(&report,param);
  int data_s=train_prepare_data(&data,data1,param);
  int t=0;
  printf("--- basic information ---\n");
  printf("ntrain,nprobe:%d,%d\n",ntrain,nprobe);
  printf("ncust,nmovie:%d,%d\n",ncust,nmovie);
  train_print_param(tmpc,param);
  printf("%s",tmpc);
  train_print_report(tmpc,report);
  printf("%s",tmpc);
  for(t=0;t<niter;t++){
    printf("--- start iteration %d ---\n",t);
    train(data,param,vector,report);
    train_print_report(tmpc,report);
    printf("%s",tmpc);
    Double score=main_judge(data2,vector,param);
    printf("SCORE: %lf\n",score);
  }
  free(tmpc);
  free(param);
  free(vector);
  free(report);
  free(data);
  free(data1);
  free(data2);
  return 0;
}
static Double main_judge(const Int64* data2, const void* vector, const void* param){
  int t;
  Double f=0;
  for(t=0;t<nprobe;t++){
    Int64 n=data2[t];
    Int n_movie=n&0xffff;
    Int n_cust=(n>>16)&0xffffff;
    Int n_ranking=(n>>40)&0xff;
    Int n_time=(n>>48)&0xffff;
    Double score=train_get_score(vector,param,n_cust,n_movie,n_time)-n_ranking;
    f+=score*score;
  }
  return sqrt(f/nprobe);
}


