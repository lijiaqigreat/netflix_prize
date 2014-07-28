
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

#define MALLOC(type,size) (type*)malloc(sizeof(type)*(size));
#define DEBUG(a,...) printf(a,##__VA_ARGS__)
/*
 * Int64:
 *  raw data:time<<48+raking<<40+cust<<16+movie<<0;
    processed:i<<32+ranking<<16+time<<0;
 */

typedef double Double;
typedef uint64_t Int64;
typedef uint32_t Int;

static const int ndata=  100480507;
static const int ncust=  480189;
static const int nmovie= 17770;
static const int nmaxid= 2649429;
static const int nprobe= 1408395;
static const int nprobem=16938;
static const int ntrain= 99072112;
static const int nrank=4;



int train_prepare_param(void** param,int argc,char** argv);
int train_prepare_vector(void **vector,void *param);
int train_prepare_report(void **report,void *param);
int train_prepare_data(void** datap,const Int64* data1,const void* param);

void train_print_param(char* rt,const void *param);
void train_print_report(char* rt,const void *report);

void train(const void* data,const void* param, void* vector,void* report_x_2);
Double train_get_score(const void* vector,const void *param,int custi,int moviei,int time);

