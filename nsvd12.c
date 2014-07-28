#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#define MALLOC(type,size) (type*)malloc(sizeof(type)*(size));
typedef double Double;
typedef uint32_t Int;
typedef struct{
  Int ndata;
  Int ncust;
  Int nmovie;
  Int ndim;
  Int *ind1;
  Int *ind2;
} Data1;
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
} Data2;

static Data1 getData1(char *path);

static Data1 getData1(char *path){
  Data1 td;
  FILE *file=fopen(path,"rb");
  printf("tmp: %d\n",(int)file);
  do{
    Int* buffer=MALLOC(Int,4);
    int n = fread(buffer, sizeof(Int), 4, file);
    td.ndata=buffer[1];
    td.ncust=buffer[2];
    td.nmovie=buffer[3];
    td.ndim=10;
    td.ind1=MALLOC(Int,td.ndata);
    td.ind2=MALLOC(Int,td.ndata);
  }while(0);
  //init custind
  fread(td.ind1,sizeof(Int)*td.ndata,1,file);
  fread(td.ind2,sizeof(Int)*td.ndata,1,file);

  fclose(file);
  return td;
}


void main_generate(char* path,int random_seed,Double percentage, void** train, void **qualify){
  srand(random_seed);
  printf("start generate\n");
  
  Data1 original = getData1(path);
  printf("read done\n");
  int qualifyn=(int)(original.ndata*percentage);
  int t1;
  for(t1=0;t1<qualifyn;t1++){
    int datai=rand()%original.ndata;
    if(original.ind1[datai] >= original.ncust){
      t1--;
      continue;
    }
    original.ind1[datai]+=original.ncust;
  }
  printf("generating data\n");
  Int *qualify_=(Int*)malloc(sizeof(Int)*(4+qualifyn*2));
  *qualify=qualify_;
  printf("test ptr: %lu\n",*qualify);
  qualify_[0]=31415;
  qualify_[1]=qualifyn;
  qualify_[2]=original.ncust;
  qualify_[3]=original.nmovie;
  qualify_+=4;
  Int* qind1=qualify_;
  Int* qind2=qualify_+qualifyn;

  int trainingn=original.ndata-qualifyn;

  printf("generating data\n");
  //init accurence count
  Int *train_=MALLOC(Int,original.ncust+original.nmovie+2+trainingn*2);
  *train=train_;
  train_[0]=trainingn;
  train_[1]=original.ncust;
  train_[2]=original.nmovie;
  train_+=3;
  Int *p1=train_;
  train_+=original.ncust+1;
  Int *p2=train_;
  train_+=original.nmovie+1;
  Int *ind1=train_;
  train_+=trainingn;
  Int *ind2=train_;
  train_+=trainingn;
  memset(p1+1,0,sizeof(Int)*original.ncust);
  memset(p2+1,0,sizeof(Int)*original.nmovie);
  qualifyn=0;
  for(t1=0;t1<original.ndata;t1++){
    //training
    if(original.ind1[t1] < original.ncust){
      p1[original.ind1[t1]+1]++;
      p2[original.ind2[t1]/5+1]++;
    }else{
      qind1[qualifyn]=original.ind1[t1]-original.ncust;
      qind2[qualifyn]=original.ind2[t1];
      qualifyn++;
    }
  }
  p1[0]=0;
  for(t1=0;t1<original.ncust;t1++){
    p1[t1+1]+=p1[t1];
  }
  p2[0]=0;
  for(t1=0;t1<original.nmovie;t1++){
    p2[t1+1]+=p2[t1];
  }
  Int *p1_=MALLOC(Int,original.ncust+1);
  memcpy(p1_,p1,sizeof(Int)*(original.ncust+1));
  Int *p2_=MALLOC(Int,original.nmovie+1);
  memcpy(p2_,p2,sizeof(Int)*(original.nmovie+1));
  for(t1=0;t1<original.ndata;t1++){
    if(original.ind1[t1]>=original.ncust){
      continue;
    }
    int id;
    id=p1_[original.ind1[t1]]++;
    ind1[id]=original.ind2[t1];
    id=p2_[original.ind2[t1]/5]++;
    ind2[id]=original.ind1[t1]*5+(original.ind2[t1]%5);
  }
  printf("DONE\n");
}
