/*******************************
Copyright (c) 2009
Univeristy of Southern California
Information Sciences Institute
Karandeep Singh
karan@east.isi.edu

Example code to show usage of
FFTW on a single tile. It also
stores generated plan in a wisdom file

To find performance of mulitple FFTWs
on simulator, use code that reads wisdom
and just executes plan in parallel
That would save time as simulation
time for creating plans can be huge.

*******************************/


#include <stdio.h>
#include <fftw3.h>
//#include <sys/profiler.h>
#include <arch/cycle.h>
#include <ilib.h>
#include <math.h>
#ifndef NUM_PROCESSES
#define NUM_PROCESSES 1
#endif

#define ILIB_GO_PARALLEL 0
#define VERIFY 1

#ifndef VERIFY
#define WARM_DCACHE 1
#define WARM_ICACHE 1
#endif

#ifdef WARM_ICACHE
#define M 2
#else
#define M 1
#endif

#ifdef VERIFY
#define N 32
#define RELATIVE_ERROR 0.01
#define ABSOLUTE_ERROR 0.01
#define INFILE_STR "32_input.dat"
#define OUTFILE_STR "32_output.dat"
#endif

#ifndef N
#define N 16
#endif

main(){
#ifdef ILIB_GO_PARALLEL
 ilib_init();
 int rank = 0;

 if (ilib_proc_go_parallel(NUM_PROCESSES, NULL) != ILIB_SUCCESS){
   printf("Failed to go_parallel().\n");
   ilib_finish();
   return -1;
  }
 // Each process has a unique rank.
 rank = ilib_group_rank(ILIB_GROUP_SIBLINGS);
 
#else
int rank = 0;
#endif
 
#ifdef VERIFY
 int correct = 1;
#endif

 int i;
 int j;
 fftw_complex *in, *out;
 fftw_plan p;
 char *wisdom_string;
 long long int start_ctr,end_ctr;
// FILE *wisdom_file;

//wisdom_file = fopen ("wisdom-file", "r");

// fftw_forget_wisdom();

// Allocate memory for input and output arrays
 in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
 out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

 
 /*
   fftw_plan fftw_plan_dft_r2c_1d(int n, double *in, fftw_complex *out,
   unsigned flags);
   fftw_plan fftw_plan_dft_c2r_1d(int n, fftw_complex *in, double *out,
   unsigned flags);
 */

// fftw_import_wisdom_from_file(wisdom_file);

//Plan for the given size
 p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_MEASURE);
 //p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_WISDOM_ONLY);

 
#ifdef VERIFY
 double *vout;
 FILE *fi, *fo;
 
 vout = (double *)malloc(sizeof(double)*2*N);
 fi = fopen (INFILE_STR, "r");
 fo = fopen (OUTFILE_STR, "r");
 
 for (i = 0; i < N; i++) {
   fscanf (fi, "%lf %lf\n", &in[i][0], &in[i][1]);
 }
 for (i = 0; i <N; i++) {
   fscanf (fo, "%lf %lf\n", &vout[2*i], &vout[2*i+1]);
 }
 
#else
 
 for(i=0;i<N;++i){
   in[i][0] = rank*i + 1.0;
   in[i][1] = 0.0;
 }
 
#endif

//warm the data caches if required 
#ifdef WARM_DCACHE
 for(i=0;i<N;++i){
   in[i][0] =5.0;
   in[i][1] = 5.0;
   out[i][0] = 0.0;
   out[i][1] = 0.0;
 }
#endif
 
#ifndef WARM_ICACHE
 start_ctr = get_cycle_count();
#endif
for(i = 0;i<N;++i)
{
out[i][0]=0.0;
out[i][1]=0.0;
}



 for(i=0; i<2; ++i){

// for warm icache, we need to execute the loop twice   
#ifdef WARM_ICACHE
   if(i==1){
     start_ctr = get_cycle_count();
   }
#endif
for (j = 0;j<N;++j){
printf("i: %d %2.15lf %2.15lf %2.15lf %2.15lf\n",i,in[j][0],in[j][1],out[j][0],out[j][1]);
}   
//Execute the plan already generated before
   fftw_execute(p); /* repeat as needed */   
 for(j = 0;j<N;++j){
printf("i: %d %2.15lf %2.15lf %2.15lf %2.15lf\n",i,in[j][0],in[j][1],out[j][0],out[j][1]); 
}
}

 end_ctr = get_cycle_count();
 fftw_print_plan(p);
 
//fftw_export_wisdom_to_file(wisdom_file);
 wisdom_string = fftw_export_wisdom_to_string();
 if(wisdom_string != NULL)
   printf("\nAccumulated wisdom: %s\n",wisdom_string);
 else
   printf("No wisdom generated!\n");
 
#ifdef VERIFY  
 for (i = 0; i < N; ++i) {
   if ((fabs(out[i][0] - vout[2*i]) > fabs(RELATIVE_ERROR*vout[2*i])) 
       && 
       (fabs(out[i][0] - vout[2*i]) > ABSOLUTE_ERROR)
       &&
       (fabs(out[i][1] - vout[2*i+1]) > fabs(RELATIVE_ERROR*vout[2*i+1])) 
       && 
       (fabs(out[i][1] - vout[2*i+1]) > ABSOLUTE_ERROR)){
     printf ("ERROR: fftw_out[%d]= %2.15lf, vout[%d]= %2.15lf\n",
	     i,out[i][0],i, vout[2*i]);
     correct = 0;
   }
 }
 if (correct) {
   printf ("PROGRAM VERIFIED: CORRECT\n");
 }
 else {
   printf ("PROGRAM NOT VERIFIED: ERROR FOUND\n");
 }
 fclose (fi);
 fclose (fo);
 
#endif

 printf("Tile#%d\tStart ctr: %lld\tEnd ctr: %lld\tCycles: %lld\n",
	rank,start_ctr,end_ctr,end_ctr-start_ctr);
 
 fftw_destroy_plan(p);
 fftw_free(in); 
 fftw_free(out);
 fftw_free(wisdom_string);
 
#ifdef ILIB_GO_PARALLEL
 ilib_finish();
#endif

}
