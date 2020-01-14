/*
 * fftw test -- double precision
 */

#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <math.h>
#include <string.h>
#include "fftw_tools.h"
 
void get_range_reference_chirp(fftw_complex *chirp,double slope, double range_pulse_length, double fs, double chirp_length, int numCols)
{
	// RETURNS CONJUGATE OF FFT OF RANGE REFERENCE CHIRP SIGNAL
	double dt = 1/fs;
	double npts = floor(range_pulse_length*fs);
	double bottom  = -npts/2;
	double pi =3.141592653589793;
	double t = 0;
	double phase = 0;
	int i = 0;
	for (i = 0; i<chirp_length;i++){
		t = bottom*dt;
		phase = pi*slope*t*t;
		chirp[i][0]=cos(phase);
		chirp[i][1]=sin(phase);
		bottom=bottom+1;
	}

	for (i = chirp_length;i<numCols ;i++){
		chirp[i][0]=0;
		chirp[i][1]=0;
		}
	fftw_plan p;
	p = fftw_plan_dft_1d(numCols,chirp,chirp,FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);
	complex_conjugate(chirp,numCols);
}

void range_compress(fftw_complex *rawdata, fftw_complex *ref, int nValid, int numRows,int numCols)
{


	int i = 0;
	int j = 0;
        
	plan_rows(rawdata,numRows,numCols,FFTW_FORWARD);

	repeat_complex_multiply(rawdata,ref,numRows,numCols);
	printf("made it past repeat multiply\n");
	plan_rows(rawdata,numRows,numCols,FFTW_BACKWARD);
	int ind = 0;
        for(i = 0;i<numRows;i++){
		for(j = 0;j<nValid;j++){
			rawdata[ind][0] = rawdata[ind][0]/numCols;
			rawdata[ind][1] = rawdata[ind][1]/numCols;
			ind = ind+1;
			}
		for(j = nValid;j<numCols;j++){
			rawdata[ind][0] = 0;
			rawdata[ind][1] = 0;
			ind = ind+1;
			}
	}
	
}

double getPhChg(fftw_complex *data, int numRows, int col, int numCols)
{

        int i = 0;
        int j = 0;
        double PhChg = 0;

        double l1_real = 0;
        double l1_imag = 0;
        double l2_real = 0;
        double l2_imag = 0;

        double l3_real = 0;
        double l3_imag = 0;
        int l1ind = 0;
        int l2ind = 0;
        double real = 0;
        double imag = 0;
        for (i = 1;i<numRows;i++){
                l1ind = i*numCols +col;
                l2ind = (i-1)*numCols + col;
                l1_real = data[l1ind][0];
                l1_imag = data[l1ind][1];
                l2_real = data[l2ind][0];
                l2_imag = -1*data[l2ind][1];
                real = l1_real*l2_real -(l1_imag*l2_imag);
                imag = l1_real*l2_imag + l1_imag*l2_real;
                l3_real = l3_real+real;
                l3_imag = l3_imag+imag;
        }
	//printf("%15.10f %15.10f\n",l3_real,l3_imag);
        PhChg = atan2(l3_imag,l3_real);
        return PhChg;


}

double get_doppler(fftw_complex *data, double prf, int numRows, int numCols, int nValid)
{
	double doppler = 0;
	int i = 0;

	for (i = 0; i<nValid;i++){
		doppler = doppler + getPhChg(data,numRows,i,numCols);
	}

	doppler = doppler/nValid;
	doppler = doppler/(2*3.1415926535);
	doppler = doppler*prf;

	return doppler;
}

void azimuth_chirp(fftw_complex *chirp, double doppler, double slope, double range_pulse_length, double fs, int numLines)
{

        double dt = 1/fs;
        double npts2 = floor(range_pulse_length*fs);
	int npts = (int) floor(range_pulse_length*fs);
        double bottom = -npts2/2;
        double pi =3.141592653589793;
        double t = 0;
        double phase = 0;
	int i = 0;
        
  
	for (i = 0;i<numLines;i++){
		chirp[i][0]=0;
		chirp[i][1]=0;
	}

 
        for (i = 0;i<npts;i++){
                t = bottom*dt;
                phase = 2*pi*doppler*t + pi*slope*t*t;
                chirp[i][0] = cos(phase);
                chirp[i][1] = sin(phase);
                bottom = bottom+1;
        }

        for (i = npts;i<numLines;i++){
                chirp[i][0]=0;
                chirp[i][1]=0;
        }


}


int get_azimuth_range(double doppler, int nValid, double fs, double ro, double vEff, double lambda, double prf,int numLines, double L)
{
	double c = 2.99792458e8;
	double range = ro + (nValid-1)*(c/(2*fs));
	double squint = sqrt(1-pow(((lambda*doppler)/(2*vEff)),2));
	double rdc = range/squint;
	double az_beamwidth = rdc*(lambda/L)*0.8;
	double tau_az = az_beamwidth/vEff;
	double nPtsAz = ceil(tau_az*prf);
	int validAzPts = 0;
	validAzPts = numLines- (int) nPtsAz;
	return validAzPts;
}

void get_azimuth_reference(fftw_complex *az_ref, double doppler, int nValid, double fs, double ro, double vEff, double lambda, double prf, int numLines, double L,int numCols)
{
	double c = 2.99792458e8;
	int i = 0;
	int j = 0;
	//double* range = new double[nValid];
	int ind = 0;
	//double *range = (double*)malloc(nValid*sizeof(double));	
	double range = 0;
	double squint = sqrt(1-pow(((lambda*doppler)/(2*vEff)),2));
	double rdc = 0;
	double fRate = 0;
	double az_beamwidth = 0;
	double tau_az = 0;
	
	fftw_complex *chirp = (fftw_complex*)malloc(numLines*sizeof(fftw_complex));

	for(i = 0;i<nValid;i++){
		range = ro+ i*c/(2*fs);
		rdc = range/squint;		
		fRate = -(2*pow(vEff,2)/lambda)/rdc;
		az_beamwidth = rdc*(lambda/L)*0.8;
		tau_az = az_beamwidth/vEff;
		azimuth_chirp(chirp,doppler,fRate,tau_az,prf,numLines);			
		for (j = 0;j<numLines;j++){
			ind = j*numCols+i;
			az_ref[ind][0] = chirp[j][0];
			az_ref[ind][1] = chirp[j][1];
		}
	}
	for(i = nValid;i<numCols;i++){
		for(j = 0;j<numLines;j++){
			ind = j*numCols+i;
			az_ref[ind][0] = 0;
			az_ref[ind][1] = 0;
		}
	}
	
	free(chirp);
	//plan_cols(az_ref,numLines,numCols,FFTW_FORWARD);
	plan_cols_2d_forward(az_ref,numLines,numCols);
        complex_conjugate(az_ref,numLines*numCols);		
}


void rcmc2(fftw_complex *data, double doppler, double lambda, double prf, double vEff, double ro, double del_sr, int validAzPts, int numLines, int nValid,int numCols)
{
	double delta_f = 0;
	double vec1 = 0;
	double vec2 = 0;
	double mid = 0;
	double offset = 0;
	int ind = 0;
	int shiftind = 0;
	int i = 0;
	int j = 0;

	for(i = 0;i<numLines;i++){
		delta_f = i*(prf/validAzPts) + doppler;
		vec1 = 1/(sqrt(1-pow((lambda*delta_f)/(2*vEff),2)))-1;
		for(j = 0;j<numCols;j++){
			vec2 = ro +j*del_sr;
			mid = vec1*vec2;
			offset = round(mid/del_sr);
			shiftind = j + (int) offset;
			if(shiftind>(nValid-1)){
				data[ind][0] = 0;
				data[ind][1] = 0;
			}
			else{
				shiftind = i*numCols + shiftind;
				data[ind][0] = data[shiftind][0];
				data[ind][1] = data[shiftind][1];
			}
			ind = ind+1;
		}
	}
}

void azimuth_compress(fftw_complex *data, fftw_complex *chirp, int numLines, int numCols)
{
	int i = 0;
       complex_multiply(data,chirp,numLines*numCols);
       //plan_cols(data,numLines,numCols,FFTW_BACKWARD);
       plan_cols_2d_backward(data,numLines,numCols); 
      //for(i = 0;i<numLines*numCols;i++){
//	data[i][0] = data[i][0]/numLines;
//	data[i][1] = data[i][1]/numLines;
//	}
}	

double sum_pixels(fftw_complex *data,int row_start,int num_rows, int col_start,int num_cols, int numCols)
{	
	int i = 0;
	int j = 0;
	double my_sum = 0;
	int ind = 0;
	for (i = row_start;i<row_start+num_rows;i++){
		for(j = col_start;j<col_start+num_cols;j++){
			ind = (i*numCols)+j;
			my_sum = my_sum + data[ind][0];
		}
	}
	return my_sum;

}

void multilook(fftw_complex *data, double **collapsed, int numLines, int nValid, int rangelooks, int azilooks,int numCols)
{

	double mysum = 0;
	printf("%2d\n",numLines);
	int rl = numLines/azilooks;
	int cl = nValid/rangelooks;
	printf("collapsed rows: %2d, collapsed cols: %2d\n",rl,cl);	
	int rowind = 0;
	int colind = 0;
	int i = 0;
	int j = 0;
	for (i = 0;i < rl;i++){
		for(j = 0;j<cl;j++){
			rowind = i*azilooks;
			colind = j*rangelooks;
			collapsed[i][j] = sum_pixels(data,rowind,azilooks,colind,rangelooks, numCols);
			mysum = mysum+collapsed[i][j];
		}
	}	 
	printf("sum at multilook output: %15.10f\n",mysum);
}


void ProcessPatch(fftw_complex *data, fftw_complex *range_chirp, fftw_complex *azimuth_chirp,double **patch, int nValid, int numRows, int numCols,double doppler, double lambda,double prf,double vEff,double ro,double del_sr,int validAzPts,int rangelooks,int azilooks)
{
	printf("starting Range Compression\n");
	range_compress(data, range_chirp, nValid, numRows,numCols);
	printf("transforming to Doppler domain\n");
	//plan_cols(data,numRows,numCols,FFTW_FORWARD);
	plan_cols_2d_forward(data,numRows,numCols);
        printf("starting Range Cell Migration Correction\n");
	rcmc2(data, doppler, lambda, prf, vEff, ro, del_sr, validAzPts, numRows, nValid, numCols);
	printf("starting azimuth compression\n");
	azimuth_compress(data,azimuth_chirp, numRows,numCols);
	printf("taking magnitude\n");
	take_magnitude(data,numRows*numCols);
	printf("starting multilook\n");
	multilook(data,patch,validAzPts,nValid,rangelooks,azilooks,numCols);
}

