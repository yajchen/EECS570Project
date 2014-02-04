#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <complex>
#include <cmath>
#include <sys/time.h>
#include "simt.h"
#include "gem5.h"
#include "time.h"

#define serial
//#define parallel

//#define OUTPUTDATA

#define x_dim 2048
#define y_dim 2048
#define stride 2

#define mean_design 128
#define std_desgin 32


typedef char uint8;
typedef unsigned short uint16;
typedef unsigned int uint32;
typedef unsigned long long uint64;

void serial_stat(uint8 **pInp,uint16 *pmean,uint16 *pstd){
  uint64 Sum=0;
  uint64 SumCube=0;
  for(int x=0;x<x_dim;x++){
    for(int y=0;y<y_dim;y++){
      Sum+=pInp[x][y];
      SumCube+=(pInp[x][y]*pInp[x][y]);
    }
  }
  
  uint32 total_pixel=x_dim*y_dim;

  uint16 tmp = (uint16)(Sum/total_pixel);
  *pmean = tmp;
  *pstd  = (uint16)sqrt((SumCube/total_pixel) - (tmp*tmp));

}

void serial_transform(uint8 **pInp,uint8 **pOut,uint16 mean,uint16 std){

float sigma_ratio = ((float)std_desgin)/((float)std);
float tmp;
for(int x=0;x<x_dim;x++){
  for(int y=0;y<y_dim;y++){
  tmp=(pInp[x][y]-mean)*sigma_ratio+mean_design;
  pOut[x][y] = (tmp<=255)? ((uint8)tmp) : 255;
  }
}

}


void column_transform(uint8 **pInp,uint8 **pOut,uint32 col,uint16 mean,uint16 std){
float sigma_ratio = ((float)std_desgin)/((float)std);
float tmp;
for(int x=0;x<x_dim;x++){
    tmp=(pInp[x][col]-mean)*sigma_ratio+mean_design;
    pOut[x][col] = (tmp<=255)? ((uint8)tmp) : 255; 
}

}


void column_stat(uint8 **pInp,uint32 tid,uint32 *pSumPixelPerThread,uint32 *pSumCubePixPerThread){
  uint32 Sum=0;
  uint32 SumCube=0; 
  for(int x=0;x<x_dim;x++){  
    for(int s=0;s<stride;s++){
     Sum+=pInp[x][tid*stride+s];
     SumCube+=(pInp[x][tid*stride+s]*pInp[x][tid*stride+s]);
    }  
  }
  *(pSumPixelPerThread+tid)  =Sum;
  *(pSumCubePixPerThread+tid)=SumCube;
}

void parallel_stat(uint8 **pInp,uint16 *pmean,uint16 *pstd){

uint32 num_threads = y_dim/stride;
uint32 *pSumPixelPerThread=new uint32[num_threads];
uint32 *pSumCubePixPerThread=new uint32[num_threads];

ar::simt_tau::par_for(num_threads, [&](size_t tid) {
     column_stat(pInp,tid,pSumPixelPerThread,pSumCubePixPerThread);
});

uint64 Sum=0;
uint64 SumCube=0;
for(uint32 tid=0;tid<num_threads;tid++){
  Sum     += pSumPixelPerThread[tid];
  SumCube += pSumCubePixPerThread[tid];
}
 
uint32 total_pixel=x_dim*y_dim;

uint16 tmp = (uint16)(Sum/total_pixel);
*pmean = tmp;
*pstd  = (uint16)sqrt((SumCube/total_pixel) - (tmp*tmp));

delete[] pSumPixelPerThread;
delete[] pSumCubePixPerThread;
}

int main(int argc, char **argv){

uint8 **pInp=new uint8*[x_dim];
uint8 **pOut=new uint8*[x_dim];
for(int x=0;x<x_dim;x++){
  *(pInp+x)=new uint8[y_dim];
  *(pOut+x)=new uint8[y_dim];
}

srand(time(NULL));
/////////generate data/////////
for(int x=0;x<x_dim;x++){
  for(int y=0;y<y_dim;y++){
    pInp[x][y] = (uint16)(rand()%256);
  }
}
fprintf(stderr,"input image size %d x %d\n",x_dim,y_dim);

#ifdef OUTPUTDATA
fprintf(stderr,"input image\n");
for(int xidx=0;xidx<x_dim;xidx++){
   for(int yidx=0;yidx<y_dim;yidx++){
    fprintf(stderr,"%d\t",pInp[xidx][yidx]);
   }
}
fprintf(stderr,"\n");
#endif

uint16 mean,std;

struct timeval tim1,tim2,tim3,tim4;
#ifdef serial
gem5::ResetStats();
//above is initial stat.txt

//start collect stat.txt for norm-stat
gem5::DumpStats();
gettimeofday(&tim1, NULL); 
serial_stat(pInp,&mean,&std);
gettimeofday(&tim2, NULL);
gem5::ResetStats();
//end collect stat.txt for norm-stat

//start collect stat.txt for norm-tran
gem5::DumpStats();
gettimeofday(&tim3, NULL); 
serial_transform(pInp,pOut,mean,std);
gettimeofday(&tim4, NULL); 
gem5::ResetStats();
//end collect stat.txt for norm-tran

gem5::DumpStats();
//below is tail stat.txt
fprintf(stderr,"serial version\n");
#endif


#ifdef parallel
gem5::ResetStats();
//above is initial stat.txt

//start collecting stat.txt for parallel norm-stat
gem5::DumpStats();
gettimeofday(&tim1, NULL);
parallel_stat(pInp,&mean,&std);
gettimeofday(&tim2, NULL);
gem5::ResetStats();
//end collecting stat.txt for parallel norm-stat

//start collecting stat.txt for parallel norm-tran
gem5::DumpStats();
gettimeofday(&tim3, NULL); 
ar::simt_tau::par_for(y_dim, [&](size_t col) {
     column_transform(pInp,pOut,col,mean,std);
});
gettimeofday(&tim4, NULL); 
gem5::ResetStats();
//end collecting stat.txt for parallel norm-tran

gem5::DumpStats();
//below is tail stat.txt
fprintf(stderr,"parallel version\n");  
#endif


double t1=tim1.tv_sec*1000000.0+tim1.tv_usec;
double t2=tim2.tv_sec*1000000.0+tim2.tv_usec; 
double t3=tim3.tv_sec*1000000.0+tim3.tv_usec;  
double t4=tim4.tv_sec*1000000.0+tim4.tv_usec;  
fprintf(stderr,"stat@@@ %f microseconds elapsed\n", t2-t1);
fprintf(stderr,"tran@@@ %f microseconds elapsed\n", t4-t3);
fprintf(stderr,"mean=%d\tstd=%d\n",mean,std);
 

fprintf(stderr,"%d\n",pOut[88][97]);


#ifdef OUTPUTDATA
fprintf(stderr,"output image\n");
for(int xidx=0;xidx<x_dim;xidx++){
   for(int yidx=0;yidx<y_dim;yidx++){
    fprintf(stderr,"%d\t",pOut[xidx][yidx]);
   }
}
fprintf(stderr,"\n");
#endif




for(int x=0;x<x_dim;x++){
  delete[] *(pInp+x);
  delete[] *(pOut+x);
}
delete[] pInp;
delete[] pOut;

return 0;
}
