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

#define T_PFA 13
#define PFA_N_TSINC_POINTS_PER_SIDE 6
#define N_RANGE 512
#define N_PULSES 512
#define PFA_NOUT_RANGE 512
#define PFA_NOUT_AZIMUTH 512

typedef struct _complex { float re, im; } complex;
/*
#define INPUT_SIZE_SMALL 1
#define INPUT_SIZE_MEDIUM 2
#define INPUT_SIZE_LARGE 3

#ifndef INPUT_SIZE
    #define INPUT_SIZE INPUT_SIZE_SMALL
#endif

#if INPUT_SIZE == INPUT_SIZE_SMALL
    #define N_RANGE (512)
    #define N_PULSES (512)
    #define BP_NPIX_X (512)
    #define BP_NPIX_Y (512)
    #define PFA_NOUT_RANGE (512)
    #define PFA_NOUT_AZIMUTH (512)
#elif INPUT_SIZE == INPUT_SIZE_MEDIUM
    #define N_RANGE (1024)
    #define N_PULSES (1024)
    #define BP_NPIX_X (1024)
    #define BP_NPIX_Y (1024)
    #define PFA_NOUT_RANGE (1024)
    #define PFA_NOUT_AZIMUTH (1024)
#elif INPUT_SIZE == INPUT_SIZE_LARGE
    #define N_RANGE (2048)
    #define N_PULSES (2048)
    #define BP_NPIX_X (2048)
    #define BP_NPIX_Y (2048)
    #define PFA_NOUT_RANGE (2048)
    #define PFA_NOUT_AZIMUTH (2048)
#else
    #error "Unhandled value for INPUT_SIZE"
#endif
*/
void read_kern1_data_file(
    complex **data,
    double *input_start_coords,
    double *input_coord_spacing,
    double *output_coords,
    float *window){
    
    FILE *fp=NULL;
    
    fp = fopen(".//input//small_kernel1_input.bin", "rb");
    size_t n;
    for(size_t p=0;p<N_PULSES;p++){
      n = fread(*(data+p), sizeof(complex), N_RANGE, fp);    
      fprintf(stderr,"%d of input complex data of row %d",n,p);             
    }

    n = fread(input_start_coords, sizeof(double), N_PULSES, fp);
    fprintf(stderr,"%d of input double of input_start_coords",n);

    n = fread(input_coord_spacing, sizeof(double), N_PULSES, fp);
    fprintf(stderr,"%d of input double of input_coord_spacing",n);      
      
    n = fread(output_coords, sizeof(double), PFA_NOUT_RANGE, fp);
    fprintf(stderr,"%d of input double of output_coords",n);      

    n = fread(window, sizeof(float), T_PFA, fp);
    fprintf(stderr,"%d of input double of window",n);      

    fclose(fp);        
  
}

int min(int a,int b){return a<b ? a:b;}
int max(int a,int b){return a>b ? a:b;}

float sinc(float x){return x==0 ? 1.0f:((float)(sin(PI*x)/(PI*x))); }

int find_nearest_range_coord(double target_coord,double input_coord_start,double input_coord_spacing,double input_coord_spacing_inv){
    return (target_coord < input_coord_start ||target_coord >= (input_coord_start + (N_RANGE-1)*input_coord_spacing)) ? -1:((int)((target_coord - input_coord_start) * input_coord_spacing_inv +0.5));
}

void sar_interp1_pixel(
    size_t idx,
    complex **resampled,
    complex **data,
    float *window,
    double *input_coords_start,
    double *input_coords_spacing,
    double *output_coords){

    
    size_t p=idx/PFA_NOUT_RANGE;
    size_t r=idx%PFA_NOUT_RANGE;

    double input_start = input_coords_start[p];
    double input_spacing = input_coords_spacing[p];
    double input_spacing_inv = 1.0 / input_spacing;
    float scale_factor = fabs(output_coords[1] - output_coords[0]) * input_spacing_inv;

    const double out_coord = output_coords[r];
    int nearest = find_nearest_range_coord(output_coords[r], input_start, input_spacing, input_spacing_inv);
    if (nearest < 0){
      resampled[p][r].re = 0.0f;
      resampled[p][r].im = 0.0f;
    }else{

            if (fabs(out_coord - (input_start + (nearest+1)*input_spacing)) < fabs(out_coord - (input_start + (nearest)*input_spacing)))
                nearest = nearest + 1;
            

            int rmin = max((nearest - PFA_N_TSINC_POINTS_PER_SIDE),0);
            int rmax = min((nearest + PFA_N_TSINC_POINTS_PER_SIDE),(N_RANGE-1));
            int window_offset = ((nearest - PFA_N_TSINC_POINTS_PER_SIDE) < 0) ? (PFA_N_TSINC_POINTS_PER_SIDE - nearest) : 0 ;

            complex accum;
            accum.re = 0.0f;
            accum.im = 0.0f;
            float sinc_arg, sinc_val, win_val;

            for (int k = rmin; k <= rmax; ++k)
            {
                win_val = window[window_offset+(k-rmin)];
                sinc_arg = (out_coord - (input_start+k*input_spacing)) * input_spacing_inv;
                sinc_val = sinc(sinc_arg);
                accum.re += sinc_val * win_val * data[p][k].re;
                accum.im += sinc_val * win_val * data[p][k].im;
            }
           resampled[p][r].re = scale_factor * accum.re;
           resampled[p][r].im = scale_factor * accum.im; 
      
    }

}

int main(int argc, char **argv){

    size_t num_resampled_elements = N_PULSES * PFA_NOUT_RANGE;

    complex **resampled = new complex*[N_PULSES];
    complex **data      = new complex*[N_PULSES];

    for(int p=0;p<N_PULSES;p++){
      *(data+p)      = new complex[N_RANGE];
      *(resampled+p) = new complex[PFA_NOUT_RANGE];
    }
    float  *window               = new float[T_PFA];
    double *input_coords_start   = new double[N_PULSES];
    double *input_coords_spacing = new double[N_PULSES];
    double *output_coords        = new double[PFA_NOUT_RANGE];
 
    read_kern1_data_file(
      data,
      input_coords_start,
      input_coords_spacing,
      output_coords,
      window);   

#ifdef serial
  for(size_t idx=0;idx<num_resampled_elements;idx++)
   sar_interp1_pixel(
    idx,
    resampled,
    data,
    window,
    input_coords_start,
    input_coords_spacing,
    output_coords);
#endif

#ifdef parallel
   
#endif


#ifdef WRITEtoFILE
  FILE *wbfp=fopen("interp1.wb","wb");
  for(size_t p=0;p<N_PULSES;p++)
    fwrite(*(resampled+p), sizeof(complex), PFA_NOUT_RANGE, wbfp);

  fclose(wbfp);
#endif
}

