#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <complex>
#include <cmath>
#include <sys/time.h>
#include "simt.h"
#define PI 3.14159265358979323846

#define serial
//#define parallel
#define WRITEtoFILE

#define WIN_SIZE 5

#define T_PFA 13
#define PFA_N_TSINC_POINTS_PER_SIDE 6
#define N_RANGE 512
#define N_PULSES 512
#define PFA_NOUT_RANGE 512
#define PFA_NOUT_AZIMUTH 512

typedef struct _complex { float re, im; } complex;
int min(int a,int b){return a<b ? a:b;}
int max(int a,int b){return a>b ? a:b;}

float sinc(float x){return x==0 ? 1.0f:((float)(sin(PI*x)/(PI*x))); }

int find_nearest_azimuth_coord(
    double target_coord,
    const double *input_coords)
{
    int left_ind, right_ind, mid_ind;
    double left_val, right_val, mid_val;

    left_ind = 0;
    right_ind = PFA_NOUT_RANGE-1;
    mid_ind = (left_ind+right_ind)/2;
    left_val = input_coords[left_ind];
    right_val = input_coords[right_ind];
    mid_val = input_coords[mid_ind];

    if (target_coord < left_val || target_coord > right_val)
    {
        return -1;
    }

    while (right_ind - left_ind > 1)
    {
        if (target_coord <= mid_val)
        {
            right_ind = mid_ind;
            right_val = mid_val;
        }
        else
        {
            left_ind = mid_ind;
            left_val = mid_val;
        }
        mid_ind = (left_ind+right_ind)/2;
        mid_val = input_coords[mid_ind];
    }

    return mid_ind;
}

void sar_interp2_col(
    size_t col,
    size_t num_sw,
    complex ** resampled,
    complex ** data,
    float *window,
    double ** input_coords,
    double *output_coords){
//    fprintf(stderr,"idx=%d\tp=%d\tr=%d\n",idx,p,r);
    float input_spacing_avg = 0.0f;
    
    for(size_t i=0;i<N_PULSES-1;++i){
      input_spacing_avg += fabs(input_coords[col][i+1] - input_coords[col][i]);
    }

    input_spacing_avg /= (N_PULSES-1);    
    float input_spacing_avg_inv = 1.0f / input_spacing_avg; 
    float scale_factor = fabs(output_coords[1] - output_coords[0]) * input_spacing_avg_inv;

    size_t p;
    double out_coord;
    int nearest,pmin,pmax,window_offset;
    float sinc_arg,sinc_val;
    complex accum;

    for(size_t sw=0;sw<num_sw;sw++){
      for(size_t wi=0;wi<WIN_SIZE;wi++){
        p=sw*WIN_SIZE+wi;
        out_coord = output_coords[p]; 
        nearest = find_nearest_azimuth_coord(out_coord, input_coords[col]);
        if(nearest<0){resampled[p][col].re=0.0f;resampled[p][col].im=0.0f;}
        else{
            if (fabs(out_coord-input_coords[col][nearest+1]) < fabs(out_coord-input_coords[col][nearest]))
            {
                nearest = nearest + 1;
            }      
            pmin = max((nearest - PFA_N_TSINC_POINTS_PER_SIDE),0);
            pmax = min((N_PULSES-1),(nearest + PFA_N_TSINC_POINTS_PER_SIDE));
            window_offset = (nearest - PFA_N_TSINC_POINTS_PER_SIDE < 0) ? (PFA_N_TSINC_POINTS_PER_SIDE - nearest):0;
            accum.re=0.0f;
            accum.im=0.0f;
 
            for (int k = pmin; k <= pmax; ++k)
            {
              sinc_arg = (out_coord - input_coords[col][k]) * input_spacing_avg_inv;
              sinc_val = sinc(sinc_arg);
              accum.re += sinc_val * window[window_offset+(k-pmin)] * data[k][col].re;
              accum.im += sinc_val * window[window_offset+(k-pmin)] * data[k][col].im; 
            }
            resampled[p][col].re = scale_factor * accum.re;
            resampled[p][col].im = scale_factor * accum.im;        
        }
      }
    }

  for(size_t p=num_sw*WIN_SIZE;p<PFA_NOUT_AZIMUTH;p++){
        out_coord = output_coords[p]; 
        nearest = find_nearest_azimuth_coord(out_coord, input_coords[col]);
        if(nearest<0){resampled[p][col].re=0.0f;resampled[p][col].im=0.0f;}
        else{
            if (fabs(out_coord-input_coords[col][nearest+1]) < fabs(out_coord-input_coords[col][nearest]))
            {
                nearest = nearest + 1;
            }      
            pmin = max((nearest - PFA_N_TSINC_POINTS_PER_SIDE),0);
            pmax = min((N_PULSES-1),(nearest + PFA_N_TSINC_POINTS_PER_SIDE));
            window_offset = (nearest - PFA_N_TSINC_POINTS_PER_SIDE < 0) ? (PFA_N_TSINC_POINTS_PER_SIDE - nearest):0;
            accum.re=0.0f;
            accum.im=0.0f;
 
            for (int k = pmin; k <= pmax; ++k)
            {
              sinc_arg = (out_coord - input_coords[col][k]) * input_spacing_avg_inv;
              sinc_val = sinc(sinc_arg);
              accum.re += sinc_val * window[window_offset+(k-pmin)] * data[k][col].re;
              accum.im += sinc_val * window[window_offset+(k-pmin)] * data[k][col].im; 
            }
            resampled[p][col].re = scale_factor * accum.re;
            resampled[p][col].im = scale_factor * accum.im;        
        }  
  } 
}


void read_kern2_data_file(
    complex **data,
    double **input_coords,
    double *output_coords,
    float *window){
    
    FILE *fp=NULL;
    fprintf(stderr,"file open\n");
                   
    fp = fopen("small_kernel2_input.bin", "rb");
    size_t n;
    for(size_t p=0;p<N_PULSES;p++){
      n = fread(*(data+p), sizeof(complex), PFA_NOUT_RANGE, fp);    
    }
    fprintf(stderr,"%d of input complex of data\n",n);             

    for(size_t p=0;p<PFA_NOUT_RANGE;p++){
      n = fread(*(input_coords+p), sizeof(double), N_PULSES, fp);
    }
    fprintf(stderr,"%d of input double of input_coords\n",n);
    
      
    n = fread(output_coords, sizeof(double), PFA_NOUT_AZIMUTH, fp);
    fprintf(stderr,"%d of input double of output_coords\n",n);      

    n = fread(window, sizeof(float), T_PFA, fp);
    fprintf(stderr,"%d of input double of window\n",n);      

    fclose(fp);        
  
}

void read_data_file(
    char *data,
    size_t num_bytes){

    size_t nread = 0;
    FILE *fp = NULL;

    fprintf(stderr,"golen file open\n");
    fp = fopen("small_golden_kernel2_output.bin", "rb");
    nread = fread(data, sizeof(char), num_bytes, fp);
    fprintf(stderr,"%d of input double of input_start_coords\n",nread);
    fclose(fp);

}

double calculate_snr(
    const complex *reference,
    complex **test){
    size_t idx;

    double num = 0.0, den = 0.0;
    for (size_t i = 0; i < PFA_NOUT_AZIMUTH; ++i)
    {
      for(size_t j=0;j<PFA_NOUT_RANGE;++j)
      {
        idx=i*PFA_NOUT_RANGE+j;
        den += (reference[idx].re - test[i][j].re) *
               (reference[idx].re - test[i][j].re);
        den += (reference[idx].im - test[i][j].im) *
               (reference[idx].im - test[i][j].im);
        num += reference[idx].re * reference[idx].re +
               reference[idx].im * reference[idx].im;
      }
    }
    fprintf(stderr,"snr: den=%f\tnum=%f\n",den,num);   
   return (den == 0) ? 140.0 : (10.0*log10(num/den));

}

int main(int argc, char **argv){
    size_t num_resampled_elements = PFA_NOUT_AZIMUTH * PFA_NOUT_RANGE;

    complex **resampled      = new complex*[PFA_NOUT_AZIMUTH];
    complex **data           = new complex*[N_PULSES];
    double  **input_coords   = new double*[PFA_NOUT_RANGE];

    for(int p=0;p<N_PULSES;p++){
      *(data+p)           = new complex[PFA_NOUT_RANGE];
    }
    for(int p=0;p<PFA_NOUT_AZIMUTH;p++){
      *(resampled+p)      = new complex[PFA_NOUT_RANGE];   
    }
    for(int p=0;p<PFA_NOUT_RANGE;p++){
      *(input_coords+p)   = new double[N_PULSES];
    }

    float  *window               = new float[T_PFA];
    double *output_coords        = new double[PFA_NOUT_AZIMUTH];
    complex *gold_resampled      = new complex[num_resampled_elements];

    size_t num_sw = PFA_NOUT_AZIMUTH/WIN_SIZE;

    read_kern2_data_file(
      data,
      input_coords,
      output_coords,
      window);   
    fprintf(stderr,"end read in input data\n");   
 
    read_data_file(
        (char *) gold_resampled,
        sizeof(complex)*num_resampled_elements);

  struct timeval tim1,tim2;
#ifdef serial
  fprintf(stderr,"serial version start\n");  
gettimeofday(&tim1, NULL);  
    
  for(size_t col=0;col<PFA_NOUT_RANGE;col++){
    sar_interp2_col(
        col,
        num_sw,
        resampled,
        data,
        window,
        input_coords,
        output_coords);     
  }

gettimeofday(&tim2, NULL);
  fprintf(stderr,"serial version end\n");  
#endif

#ifdef parallel
  fprintf(stderr,"parallel version start\n");  
gettimeofday(&tim1, NULL);
    ar::simt_tau::par_for(PFA_NOUT_RANGE, [&](size_t col) {
      sar_interp2_col(col,num_sw,resampled,data,window,input_coords,output_coords); 
    });   
gettimeofday(&tim2, NULL);
  fprintf(stderr,"parallel version end\n");     
#endif

double t1=tim1.tv_sec+(tim1.tv_usec/1000000.0);
double t2=tim2.tv_sec+(tim2.tv_usec/1000000.0); 

fprintf(stderr,"interp1@@@ %f seconds elapsed\n", t2-t1);

#ifdef WRITEtoFILE
  FILE *wbfp=fopen("interp2.wb","wb");
  for(size_t p=0;p<PFA_NOUT_AZIMUTH;p++)
    fwrite(*(resampled+p), sizeof(complex), PFA_NOUT_RANGE, wbfp);

  fclose(wbfp);
#endif  


  double snr = calculate_snr((complex *) gold_resampled,resampled);
  fprintf(stderr,"snr=%f\n",snr);     


  for(int p=0;p<N_PULSES;p++){
    delete[] *(data+p);
  }

  for(int p=0;p<PFA_NOUT_AZIMUTH;p++){
    delete[] *(resampled+p);
  }

  for(int p=0;p<PFA_NOUT_RANGE;p++){
    delete[] *(input_coords+p);
  }

  delete[] data;
  delete[] resampled;
  delete[] window;
  delete[] input_coords;
  delete[] output_coords;
  delete[] gold_resampled;

}
