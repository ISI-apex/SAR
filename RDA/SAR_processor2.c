/*
 * fftw test -- double precision
 */

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
   int numLines = 2048;
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
   //int **arr = (int **)malloc(r*sizeof(int*));
   //for (i = 0;i<r;i++)
     //   arr[i]=(int*)malloc(col*sizeof(int));


  // int ind = 0;
  // for(i = 0;i<r;i++){
    //    for(j = 0;j<col;j++){
     //   arr[i][j]=buffer[ind];
       // ind = ind+1;
   // }
  // }
   fclose (pFile);



fftw_complex *csar = (fftw_complex*)malloc(numLines*numCols*sizeof(fftw_complex));

build_array_for_processing2(buffer,csar,2,numLines,numCols,bytesperline,header);




//   fftw_complex **csar = (fftw_complex **)malloc(lSize*sizeof(fftw_complex*));
  //  for (i = 0;i<r;i++)
    //     csar[i]=(fftw_complex*)malloc(numCols*sizeof(fftw_complex));

 //   int colind =0;
 //   for (i = 0;i<r;i++){
   //      colind = 0;
     //     for(j = numHdr;j<col;j++){
//		  if(j%2==0){
  //                csar[i][colind][0]=arr[i][j];
//		  }
  //                else{
    //              csar[i][colind][1]=arr[i][j];
      //            colind = colind+1;
        //          }
         //  }
  //  }

  //  fftw_complex* csar2 = new fftw_complex[numLines*numCols];
 //  fftw_complex *csar2 = (fftw_complex*)malloc(numLines*numCols*sizeof(fftw_complex));



  //  build_array_for_processing(csar,csar2,2,numLines,numCols);

   // fftw_complex* range_chirp = new fftw_complex[nValid];
    fftw_complex *range_chirp = (fftw_complex*)malloc(nValid*sizeof(fftw_complex));
   get_range_reference_chirp(range_chirp,slope,tau,fs,chirp_length,nValid);

  //  fftw_complex* compressed = new fftw_complex[numLines*nValid];
   fftw_complex *compressed = (fftw_complex*)malloc(numLines*nValid*sizeof(fftw_complex));
 
   range_compress(csar,range_chirp,compressed,nValid,numLines,numCols);

   double doppler = get_doppler(compressed,prf,numLines,nValid);
   int validAzPts = get_azimuth_range(doppler, nValid,fs,ro,vEff,lambda,prf,numLines,L);
   printf("doppler: %15.10f\n",doppler);
   free(csar);
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
	   
	//fftw_complex  **all_slc = (fftw_complex **)malloc(patches*validAzPts*(sizeof(fftw_complex*)));
	//for(int i = 0;i<patches*validAzPts;i++)
          //   all_slc[i]=(fftw_complex*)malloc(nValid*sizeof(fftw_complex));

printf("built all_slc\n");

   // fftw_complex* compressed_FT = new fftw_complex[numLines*nValid];
    fftw_complex *compressed_FT = (fftw_complex*)malloc(numLines*nValid*sizeof(fftw_complex));  
    
    //fftw_complex* azi_chirp = new fftw_complex[numLines*nValid];
    fftw_complex *azi_chirp = (fftw_complex*)malloc(numLines*nValid*sizeof(fftw_complex));
  
  //fftw_complex* rcmc_data = new fftw_complex[nValid*numLines];
  fftw_complex *rcmc_data = (fftw_complex*)malloc(nValid*numLines*sizeof(fftw_complex));

//    fftw_complex* slc = new fftw_complex[numLines*nValid];
  fftw_complex *slc = (fftw_complex*)malloc(numLines*nValid*sizeof(fftw_complex));  


   int valid = get_azimuth_reference(azi_chirp,doppler,nValid,fs,ro,vEff,lambda,prf,numLines,L);
	printf("%2d\n",valid);
l =0;
int slc_start = 0;
int slcind = 0;


int azilook = 20;
int rangelook = 4;

int finalrow = validAzPts/azilook;
int fr = validAzPts/azilook;
finalrow = finalrow*patches;
int finalcol = nValid/rangelook;

double **final_image = (double **)malloc(finalrow*(sizeof(double*)));
for(i = 0;i<finalrow;i++)
	final_image[i] = (double*)malloc(finalcol*sizeof(double));

fftw_complex **myslc =(fftw_complex **)malloc(validAzPts*(sizeof(fftw_complex*)));
for(i = 0;i<validAzPts;i++)
	myslc[i] = (fftw_complex*)malloc(nValid*sizeof(fftw_complex));

double **magnitude = (double **)malloc(validAzPts*(sizeof(double*)));
for(i = 0;i<validAzPts;i++)
	magnitude[i]=(double*)malloc(nValid*(sizeof(double)));

double **collapsed = (double **)malloc(finalrow*(sizeof(double*)));
for(i = 0;i<finalrow;i++)
	collapsed[i] = (double*)malloc(finalcol*(sizeof(double)));


fftw_complex **csar_all = (fftw_complex **)malloc(patches*sizeof(fftw_complex*));
for(i = 0;i<patches;i++)
	csar_all[i] = (fftw_complex*)malloc(numCols*numLines*sizeof(fftw_complex));

startind = 0;
for(i = 0;i<patches;i++){
	startind = i*validAzPts +2;
	build_array_for_processing2(buffer,csar_all[i],startind,numLines,numCols,bytesperline,header);
	//brow = brow+numLines;
}





printf("number of rows: %2d, number of cols:%2d\n",finalrow,finalcol);

int final_ind = 0;
double checksum = 0;
clock_t begin;
clock_t end;
double t;

for (i = 0; i < patches;i++)
{
	printf("patch %2d\n",i);
	checksum =0;
	slcind = 0;
	startind = l*validAzPts+2;
	printf("building array\n");
	//begin = clock();
	//build_array_for_processing(csar,csar2,startind,numLines,numCols);
        	
       //end = clock();
	//t = (end-begin)/CLOCKS_PER_SEC;
	//printf("array build time: %15.10f s\n",t); 
	printf("range compression\n");
	//begin = clock();
	range_compress(csar_all[i],range_chirp,compressed,nValid,numLines,numCols);
	//end = clock();
	//t = (end-begin)/CLOCKS_PER_SEC;
	//printf("range compression time: %15.10f s\n",t);
	printf("convert to doppler domain\n");
	//begin = clock();
	plan_cols(compressed,compressed_FT,numLines,nValid,FFTW_FORWARD);
	//end = clock();
	//t = (end-begin)/CLOCKS_PER_SEC;
	//printf("conversion time: %15.10f s\n",t);
	printf("range cell migration correction\n");
	//begin = clock();
	rcmc2(compressed_FT,rcmc_data,doppler,lambda,prf,vEff,ro,del_sr,validAzPts,numLines,nValid);
	//end = clock();
	//t = (end-begin)/CLOCKS_PER_SEC;
	//printf("range cell migration correction time: %15.10f s\n",t);
	printf("azimuth compression\n");
	//begin = clock();
	azimuth_compress(rcmc_data,azi_chirp,slc,numLines,nValid);
	//end = clock();
	//t = (end-begin)/CLOCKS_PER_SEC;
	//printf("azimuth compression time: %15.10f s\n",t);

	for (j = 0;j<validAzPts;j++){
		for(k = 0;k<nValid;k++){
			myslc[j][k][0]=slc[slcind][0];
			myslc[j][k][1]=slc[slcind][1];
			slcind = slcind+1;
		}
	}
	printf("taking magnitude\n");
	//begin = clock();
	take_magnitude(myslc, magnitude, validAzPts, nValid);
	//end = clock();
	//t = (end-begin)/CLOCKS_PER_SEC;
	//printf("magnitude time: %15.10f s\n",t);
	printf("multilook\n");
	//begin = clock();
	multilook(magnitude,collapsed,validAzPts,nValid,rangelook,azilook);
	//end = clock();
	//t = (end-begin)/CLOCKS_PER_SEC;
	//printf("multilook time: %15.10f\n",t);
	printf("saving patch\n");
	//begin = clock();
	for(j = 0 ; j < fr;j++){
		for(k = 0;k<finalcol;k++){
			final_image[final_ind][k]=collapsed[j][k];
			checksum = checksum+collapsed[j][k];
		}
		final_ind = final_ind+1;
	} 
	//end = clock();
	//t = (end-begin)/CLOCKS_PER_SEC;
	//printf("patch save time: %15.10f\n",t);
	l=l+1;
}


double rs = 0;
double is = 0;

for (i = 0;i<finalrow;i++){
	for (j = 0;j<finalcol;j++){
		rs = final_image[i][j]+rs;
	}
}

printf("sum of final image: %15.10f\n",rs);

FILE *ofile = fopen("SAR_out.bin","wb");
for(i = 0;i<finalrow;++i)
	fwrite(final_image[i],sizeof(final_image[i][0]),finalcol,ofile);

fclose(ofile);

 
//    plan_cols(compressed,compressed_FT,numLines,nValid,FFTW_FORWARD);


    //printf("%15.10f\n",doppler);
//    fftw_complex* azi_chirp = new fftw_complex[numLines*nValid];
    //printf("getting azi chirp\n");
   // int valid = get_azimuth_reference(azi_chirp,doppler,nValid,fs,ro,vEff,lambda,prf,numLines,L);
    //int valid2 = get_azimuth_range(doppler, nValid, fs, ro, vEff, lambda,prf,numLines,L);
      

 //printf("%2d %2d\n",validAzPts,valid2);
   // printf("starting rcmc\n");
  //  fftw_complex* rcmc_data = new fftw_complex[nValid*numLines];
  //  rcmc(compressed_FT,rcmc_data,doppler,lambda,prf,vEff,ro,del_sr,validAzPts,numLines,nValid);
//int ind = 0;
//    for(int i = 0;i<numLines*nValid;i++){
//	printf("%2d %15.10f %15.10f\n",i,rcmc_data[i][0],rcmc_data[i][1]);
  //  }


//    fftw_complex* slc = new fftw_complex[numLines*nValid];
  //  printf("azimuth compression\n");
   // azimuth_compress(rcmc_data,azi_chirp,slc,numLines,nValid);
    

  //  for(int i = 0;i <validAzPts*nValid;i++)
   // {
   //  printf("%2d %15.10f %15.10f\n",i,slc[i][0],slc[i][1]);
   // }
		
    
//free(csar);
//free(csar2);
//delete csar2; 
//delete range_chirp;
free(range_chirp);
//delete compressed;
free(compressed);
//delete azi_chirp;
free(azi_chirp);
//delete rcmc_data;
free(rcmc_data);
//delete slc;
free(slc);
free(magnitude);
free(final_image);
free(collapsed);
free(buffer);
//free(arr);
free(compressed_FT);
free(csar_all); 



    return 0;
}
