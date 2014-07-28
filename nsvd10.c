#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

typedef double Double;
typedef uint32_t Int;

void main_train(int iter,void* dat,void *vector);
void main_generate(char* path,int random_seed,Double percentage, void** train, void **qualify);
void main_judge(void* qualify,void* vector);
//main path
//seed
//percentage
//iter
static void* qualify;

int main(int argc,char** argv){
  if(argc!=5){
    printf("usage: ./main.out\n path of main.dat\n random seed\n percentage of qualify\n number of training iteration\n");
  }
  char* path=argv[1];
  int random_seed=atoi(argv[2]);
  Double percentage=(Double)atoi(argv[3])/100;
  int niter=atoi(argv[4]);
  
  
  void* train;
  main_generate(path,random_seed,percentage,&train,&qualify);
  main_train(atoi(argv[4]),train,0);
  
}
void request_judge(void* vector){
  main_judge(qualify,vector);
}

