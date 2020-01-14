/*
 * fftw test -- double precision
 */
#include <stdlib.h>
#include <stdio.h>
#include <fftw3.h>
#include <math.h>
#include "fftw_tools.h"
#include "SAR_tools.h"
#define N 8

int main(int argc, char *argv[])
{
   // double in1[] = { 0.00000, 0.12467, 0.24740, 0.36627,
    //                 0.47943, 0.58510, 0.68164, 0.76754
    //};
    //fftw_complex* data = new fftw_complex[2*4];
    //fftw_complex* out = new fftw_complex[2*4]; 
  // fftw_complex *data = (fftw_complex*)malloc(8*sizeof(fftw_complex));  
   fftw_complex data[8];
  int rows = 2;
    int cols = 4;
     
    //p1 = fftw_plan_many_dft(1,&N2,howmany,in,NULL,istride,idist,out,NULL,ostride,odist,FFTW_FORWARD,FFTW_MEASURE);  
   
    data[0][0]=0;
    data[1][0]=1;
    data[2][0]=2;
    data[3][0]=3;
    data[4][0]=4;
    data[5][0]=5;
    data[6][0]=6;
    data[7][0]=7;
    
    data[0][1] =7;
    data[1][1]=6;
    data[2][1]=5;
    data[3][1]=4;
    data[4][1]=3;
    data[5][1]=2;
    data[6][1]=1;
    data[7][1]=0;

int i =  0;   
    for (i = 0;i <N;i++){
    printf("%2d %15.10f %15.10f \n",i,data[i][0],data[i][1]);
    }
    printf("take the transform!\n");
   plan_rows(data,rows,cols,FFTW_FORWARD);
    //double in2[N];
    //int howmany = 2;  
    //fftw_complex  out[N];
    //fftw_plan     p1, p2;
   //int N2 = 4;
    //p1 = fftw_plan_dft_r2c_1d(N, in1, out, FFTW_ESTIMATE);
    //p1 = fftw_plan_many_dft(1,&N2,howmany,in,NULL,1,N2,out,NULL,1,4,FFTW_FORWARD,FFTW_MEASURE); 
    //p2 = fftw_plan_dft_c2r_1d(N, out, in2, FFTW_ESTIMATE);

    //fftw_execute(p1);
    for (i = 0; i < N;i++){
    printf("%2d %15.10f %15.10f \n",i,data[i][0],data[i][1]);
    }
   printf("take the inverse transform!\n");    
    plan_rows(data,rows,cols,FFTW_BACKWARD);
for (i = 0;i<8;i++){
   data[i][0] = data[i][0]/4;
   data[i][1] = data[i][1]/4;
   }   


for (i = 0; i < N;i++){
     printf("%2d %15.10f %15.10f \n",i,data[i][0],data[i][1]);
    }  

printf("take the forward column transform!\n");
    plan_cols(data,rows,cols,FFTW_FORWARD);

      for (i = 0; i < N;i++){
    printf("%2d %15.10f %15.10f \n",i,data[i][0],data[i][1]);
    }

  plan_cols(data,rows,cols,FFTW_BACKWARD);

  for (i = 0;i<8;i++){
   data[i][0] = data[i][0]/2;
   data[i][1] = data[i][1]/2;
   } 
printf("take the backward column transform!\n");
  for (i = 0; i < N;i++){
    printf("%2d %15.10f %15.10f \n",i,data[i][0],data[i][1]);
    }

printf("now calculate it with other function\n");
  plan_cols_2d_forward(data,rows,cols);


  for (i = 0; i < N;i++){
    printf("%2d %15.10f %15.10f \n",i,data[i][0],data[i][1]);
    }
printf("and take it backwards one more time\n");
plan_cols_2d_backward(data,rows,cols);
 for (i = 0; i < N;i++){
    printf("%2d %15.10f %15.10f \n",i,data[i][0],data[i][1]);
    }
  
//    fftw_complex* ta1 = new fftw_complex[2];
  //  fftw_complex* ta2 = new fftw_complex[2];

 //  ta1[0][0]=1;
 //  ta1[0][1]=2;
 //  ta1[1][0]=3;
//   ta1[1][1]=4;
//   ta2[0][0]=1;
 //  ta2[0][1]=-2;
  //ta2[1][0]=3;
  // ta2[1][1]=-4;
   //   printf("test multiply!\n");
 //  complex_multiply(ta1,ta2,2);
 // for (int i = 0;i<2;i++){
   //  printf("%2d %15.10f %15.10f \n",i,ta1[i][0],ta1[i][1]);
//}

  // double tau = 3.712e-5;
 //  double fs = 1.89625e7;
 //  double slope = 4.19e11;
 //  double numBytes = 11644;
  // double numHdr = 412;
  // int nValid = int((numBytes-numHdr)/2 - round(tau*fs));
  // double chirp_length = floor(tau*fs);

   //printf("%15.10f %15.10f %15.10f %15.10f
   
   
//   fftw_complex* chirp = new fftw_complex[nValid];
    
 //  get_range_reference_chirp(chirp,slope,tau,fs,chirp_length,nValid);


  // for (int i = 0;i<nValid;i++){
   // printf("%2d %15.10f %15.10f \n",i,chirp[i][0],chirp[i][1]);
    //}

   




    return 0;
}
