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
   
 fftw_complex *data = (fftw_complex*)malloc(22000*5616*sizeof(fftw_complex));  

  int rows = 2048;
    int cols = 5616;
     
    //p1 = fftw_plan_many_dft(1,&N2,howmany,in,NULL,istride,idist,out,NULL,ostride,odist,FFTW_FORWARD,FFTW_MEASURE);  
   
int i = 0;
for(i = 0;i<rows*cols;i++){
data[i][0] = i;
data[i][1] = i;
}




   my_transpose(data,cols,rows);
   plan_rows(data,cols,rows,FFTW_FORWARD);
   plan_rows(data,cols,rows,FFTW_BACKWARD);

}
