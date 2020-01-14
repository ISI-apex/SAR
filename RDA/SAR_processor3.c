



#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <math.h>
#include <time.h>
#include "fftw_tools.h"
#include "SAR_tools.h"
#define N 8



int main(int argc, char *argv[])
{



   double tau = 3.712e-5;
   double fs = 1.89625e7;
   double slope = 4.19e11;
   double numBytes = 11644;
   int bytesperline = 11644;
   double numHdr = 412;
   int header = 412;
   int numLines = 20000;
   int numCols = 5616;
   double c = 2.99792458e8;
   int nValid = (int) ((numBytes-numHdr)/2 - round(tau*fs));
   double L = 10;
   double v = 7125;
   double Re = 6378144;
   double h = 790000;
   double vEff = v*sqrt(Re/(Re+h));
   double ro = 8.281220820121720e5;
   double lambda = 0.056666;
   double prf = 1679.902394;
   double del_sr = c/(2*fs);
   double chirp_length = floor(tau*fs);
   int i = 0;
   int j = 0;
   int k = 0;

  FILE * pFile;
  long lSize;
  char * buffer;
  size_t result;

  pFile = fopen ( "myfile.bin" , "rb" );
  if (pFile==NULL) {fputs ("File error",stderr); exit (1);}

  // obtain file size:
  fseek (pFile , 0 , SEEK_END);
  lSize = ftell (pFile);
  rewind (pFile);
 // printf("%d\n",lSize/11644);
  // allocate memory to contain the whole file:
  buffer = (char*) malloc (sizeof(char)*lSize);
  if (buffer == NULL) {fputs ("Memory error",stderr); exit (2);}

  // copy the file into the buffer:
  result = fread (buffer,1,lSize,pFile);
  if (result != lSize) {fputs ("Reading error",stderr); exit (3);}

  /* the whole file is now loaded in the memory buffer. */

  //printf("%d\n",buffer[1000*1500-1]);
        // terminate

   int r = lSize/11644;
   int col = 11644;
   printf("r=%2d\n",r); 
   fclose (pFile);



//fftw_complex *csar = (fftw_complex*)malloc(numLines*numCols*sizeof(fftw_complex));
fftw_complex *csar;

csar = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*numLines*numCols);



printf("building array for processing\n");
build_array_for_processing2(buffer,csar,2,numLines,numCols,bytesperline,header);
//fftw_complex *range_chirp = (fftw_complex*)malloc(numCols*sizeof(fftw_complex));
fftw_complex *range_chirp;
range_chirp = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*numCols);
printf("building range reference chirp\n");
get_range_reference_chirp(range_chirp,slope,tau,fs,chirp_length,numCols);
printf("range compression\n");
range_compress(csar,range_chirp,nValid,numLines,numCols);
printf("getting doppler\n");
double doppler = get_doppler(csar,prf,numLines,numCols,nValid);
int validAzPts = get_azimuth_range(doppler, nValid,fs,ro,vEff,lambda,prf,numLines,L);
printf("doppler: %15.10f\n",doppler);
fftw_free(csar);
printf("validAzPts:%2d\n",validAzPts);
int patches = 0;
int endind = 0;
int l =0;
int startind = 0;

while(endind<r)
{
	startind = l*validAzPts+2;
	endind = startind+(numLines-1);
	l = l+1;
	patches = patches+1;
}
	patches = patches-1;
	printf("patches=%2d\n",patches);
	   


//fftw_complex *azi_chirp = (fftw_complex*)malloc(numLines*numCols*sizeof(fftw_complex));
fftw_complex *azi_chirp;
azi_chirp = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*numLines*numCols);


get_azimuth_reference(azi_chirp,doppler,nValid,fs,ro,vEff,lambda,prf,numLines,L,numCols);

//for(i = 0;i<numLines*numCols;i++){
//printf("%2d %15.10f %15.10f\n",i,azi_chirp[i][0],azi_chirp[i][1]);
//}


int azilooks = 20;
int rangelooks = 4;

int finalrow = validAzPts/azilooks;
int finalcol = nValid/rangelooks;

double ***final_image = (double ***)malloc(patches*(sizeof(double**)));
for(i = 0;i<patches;i++){
	final_image[i] = (double **)malloc(finalrow*(sizeof(double*)));
	for (j = 0;j<finalrow;j++){
		final_image[i][j] = (double*)malloc(finalcol*(sizeof(double)));
	}
}
			
fftw_complex **csar_all = (fftw_complex **)malloc(patches*sizeof(fftw_complex*));
for(i = 0;i<patches;i++)
	csar_all[i] = (fftw_complex*)malloc(numCols*numLines*sizeof(fftw_complex));

startind = 0;
for(i = 0;i<patches;i++){
	startind = i*validAzPts +2;
	build_array_for_processing2(buffer,csar_all[i],startind,numLines,numCols,bytesperline,header);
}





printf("number of rows: %2d, number of cols:%2d\n",finalrow,finalcol);


for (i = 0; i < patches;i++)
{
	ProcessPatch(csar_all[i], range_chirp,azi_chirp,final_image[i], nValid, numLines, numCols, doppler, lambda, prf, vEff, ro, del_sr, validAzPts, rangelooks, azilooks);
}

printf("finished processing patches\n");
double rs = 0;
double is = 0;

for (i = 0;i<patches;i++){

	for (j = 0;j<finalrow;j++){
		for (k = 0;k<finalcol;k++){
		rs = final_image[i][j][k]+rs;
		}
	}
	
}

printf("sum of final image: %15.10f\n",rs);

FILE *ofile = fopen("SAR_out.bin","wb");
for(i = 0;i<patches;++i){
	for(j = 0;j<finalrow;++j){
		fwrite(final_image[i][j],sizeof(final_image[i][0]),finalcol,ofile);
	}
}

fclose(ofile);
fftw_free(range_chirp);
fftw_free(azi_chirp);
free(final_image);
free(buffer);
fftw_free(csar_all); 
return 0;
}
