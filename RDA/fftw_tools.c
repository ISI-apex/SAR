/*
 * fftw test -- double precision
 */

#include <stdio.h>
#include <fftw3.h>
#include <math.h>

void plan_rows(fftw_complex *data,int rows, int cols, int direction)
{
    fftw_plan p1;
    int howmany = rows;
    int N2 = cols;
    int istride = 1;
    int ostride = 1;
    int idist = cols;
    int odist = cols;
     
    p1 = fftw_plan_many_dft(1,&N2,howmany,data,NULL,istride,idist,data,NULL,ostride,odist,direction,FFTW_ESTIMATE);   
    fftw_execute(p1);

    fftw_destroy_plan(p1);
}

void plan_cols(fftw_complex *data, int rows, int cols, int direction)
{
    fftw_plan p1;
    int howmany = cols;
    int N2 = rows;
    int idist = 1;
    int odist = 1;
    int istride = cols;
    int ostride = cols;

   p1 = fftw_plan_many_dft(1,&N2,howmany,data,NULL,istride,idist,data,NULL,ostride,odist,direction,FFTW_ESTIMATE);
   fftw_execute(p1);
   fftw_destroy_plan(p1);
}

void plan_cols_2d_forward(fftw_complex *data,int rows, int cols)
{
   fftw_plan p1;
   p1 = fftw_plan_dft_2d(rows,cols,data,data,FFTW_FORWARD,FFTW_ESTIMATE);
   fftw_execute(p1);
   fftw_destroy_plan(p1);
   plan_rows(data,rows,cols,FFTW_BACKWARD);
   int i = 0;
   for(i=0;i<rows*cols;i++)
   {
    data[i][0]=data[i][0]/cols;
    data[i][1]=data[i][1]/cols;
   }


}
void plan_cols_2d_backward(fftw_complex *data,int rows, int cols)
{
  fftw_plan p1;
  p1 = fftw_plan_dft_2d(rows,cols,data,data,FFTW_BACKWARD,FFTW_ESTIMATE);
 fftw_execute(p1);
 fftw_destroy_plan(p1);
 plan_rows(data,rows,cols,FFTW_FORWARD);
 int i = 0;
 for (i =0;i<rows*cols;i++)
{
   data[i][0]=(data[i][0]/(rows*cols));
  data[i][1] = (data[i][1]/(rows*cols));
}
}


void complex_multiply(fftw_complex *array1, fftw_complex *array2,int size)
{

double real = 0;
double imag = 0;
int i  =0;
      for (i =0; i<size;i++){
      real = array1[i][0]*array2[i][0] - array1[i][1]*array2[i][1];
      imag = array1[i][0]*array2[i][1] + array1[i][1]*array2[i][0];
      array1[i][0] = real;
      array1[i][1] = imag;
      }

}

void complex_conjugate(fftw_complex *array1,int size)
{
    int i = 0;
	for (i = 0;i<size;i++){
      array1[i][1]=-1*array1[i][1];
    }
}


void repeat_array(fftw_complex *array1, fftw_complex *array2,int numRows, int numCols)
{
	int ind = 0;
	int i = 0;
	int j = 0;
	for (i = 0; i < numRows ;i++){
		for (j = 0; j < numCols;j++){
			array1[ind][0]=array2[j][0];
			array1[ind][1]=array2[j][1];
			ind = ind+1;
			}
	}
}

void build_array_for_processing(fftw_complex **fullarray, fftw_complex *subarray, int start_row,int numLines,int numCols)
{
	double real = 0;
	double imag = 0;
	int ind = 0;
	int i = 0;
	int j = 0;
	for(i = start_row ;i< (start_row+numLines);i++){
		for (j = 0;j<numCols;j++){
			subarray[ind][0]=fullarray[i][j][0];
			subarray[ind][1]=fullarray[i][j][1];
			real = real+fullarray[i][j][0];
			imag = imag+fullarray[i][j][1];
			ind = ind+1;
		}
	}
	real = real/(numLines*numCols);
	imag = imag/(numLines*numCols);

	for (i = 0; i<(numLines*numCols);i++){
		subarray[i][0]=subarray[i][0]-real;
		subarray[i][1]=subarray[i][1]-imag;
	}

}

void build_array_for_processing2(char *filebuffer, fftw_complex *subarray,int start_row, int numLines,int numCols,int bytesperline, int header)
{
double real = 0;
double imag = 0;
int ind = 0;
int find = 0;
int i = 0;
int j = 0;

for(i = start_row;i<(start_row+numLines);i++){
	for (j = header;j<(2*numCols+header);j++){
		find = (i*bytesperline)+j;
		if(find%2==0){
		subarray[ind][0]=filebuffer[find];
		real = real+subarray[ind][0];
		}
		else{
		subarray[ind][1]=filebuffer[find];
		imag = imag+subarray[ind][1];
		ind = ind+1;
		}
	}
}
real = real/(numLines*numCols);
imag = imag/(numLines*numCols);

for(i = 0;i<(numLines*numCols);i++)
{
	subarray[i][0]=subarray[i][0]-real;
	subarray[i][1]=subarray[i][1]-imag;
}


}




void take_magnitude(fftw_complex *complex_array, int size)
{
	double real = 0;
	double imag = 0;
	int i = 0;
	for (i = 0;i<size;i++){
			real = complex_array[i][0];
			imag = complex_array[i][1];
			complex_array[i][0] = sqrt(pow(real,2)+pow(imag,2));
	}
}


void repeat_complex_multiply(fftw_complex *long_array,fftw_complex *short_array,int numRows,int numCols)

{
double real = 0;
double imag = 0;
int i = 0;
int j = 0;
int ind = 0;

for (i = 0;i<numRows;i++){
	for(j = 0;j<numCols;j++){
		real = long_array[ind][0]*short_array[j][0] - long_array[ind][1]*short_array[j][1];
		imag = long_array[ind][0]*short_array[j][1] + long_array[ind][1]*short_array[j][0];
		long_array[ind][0] = real;
		long_array[ind][1] = imag;
		ind = ind+1;
	}
}




}


void my_transpose(fftw_complex *data, int cols, int rows)

{
	int start, next, i;
	double real_tmp;
	double imag_tmp;

	
	for(start = 0;start<=cols*rows-1;start++){
		next = start;
		i = 0;
		do{ i++;
			next = (next%rows)*cols + next/rows;
		} while (next>start);
		if(next<start || i ==1) continue;
		
		real_tmp = data[next = start][0];
		imag_tmp = data[next = start][1];
		do {
			i = (next%rows)*cols+next/rows;
			data[next][0] = (i== start) ? real_tmp : data[i][0];
			data[next][1] = (i ==start)? imag_tmp : data[i][1];
		} while (next > start);
	}
}
	





