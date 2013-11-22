//
// reference Adrea Vedaldi
//
#include "mex.h"
#include <stdlib.h>
#include<math.h>
#include <string.h>

#define BINS_SIZE 36
#define PI 3.14159265358979

#define min(a,b)     (((a)<(b))?(a):(b))
#define max(a,b)     (((a)>(b))?(a):(b))

void mexFunction(int nout,mxArray *out[],int nin,const mxArray *in[])
{
	int S0 = 3;
	double sigma0 = 1.6;
	const double factor = 1.5; 
	double Hist_pt[BINS_SIZE];
	enum {IN_LOC=0,IN_DOGS,IN_S,IN_SIGMA} ;
	enum {OUT_Q=0};
	if(nin<2)
	{
		mexErrMsgTxt("At Least two input argument is required.");
	}
	if (nin>=3)
	{
		S0 = *mxGetPr(in[IN_S]) ;		
	}
	if (nin >=4)
	{
		sigma0 = *mxGetPr(in[IN_SIGMA]);
	}
	// DoGs 
	int M,N,S;
	const int* dims = mxGetDimensions(in[IN_DOGS]);
	M = dims[0] ; //y range [1,M]
	N = dims[1] ; //x range [1,N]
	S = dims[2] ; //scale   [1,S]
	if(S < 3 || M < 3 || N < 3) {
		mexErrMsgTxt("All dimensions of DOG must be not less than 3.") ;
	}
	const double* DoGs_ptr = mxGetPr(in[IN_DOGS]) ; // head pointer of DoGs vectors (M*N*S)
	// Local Max Points	
	const double* Loc_Max_ptr = mxGetPr(in[IN_LOC]) ; // head pointer of Local Max vector (3*LOCAL_MAX_NUM)
	int LOCAL_MAX_NUM = mxGetN(in[IN_LOC]) ; //the number of Local Maximum points .
	if( LOCAL_MAX_NUM == 0) {
		out[OUT_Q] = mxCreateDoubleMatrix(4, 0, mxREAL) ;
		return ;
	}

////////////////////////////////////////////////////////////
	// (M,N,S)
	const int Y_OFFSET = 1;
	const int X_OFFSET = M;
	const int S_OFFSET = M*N ;

	int buf_size = LOCAL_MAX_NUM*(3+1);
	double* buf_head = (double*) mxMalloc(buf_size*sizeof(double)) ;
	double* buf_it = buf_head;
	double* buf_tail = buf_head + buf_size;

	// for every key point 
	for (int p = 0;p<LOCAL_MAX_NUM;p++)
	{
		const double x = *Loc_Max_ptr++;
		const double y = *Loc_Max_ptr++;
		const double s = *Loc_Max_ptr++;
		int xi = (int)(x+0.5) ;
		int yi = (int)(y+0.5);		
		int si = (int)(s+0.5);
		if (xi<0 || xi>N-1 ||
			yi<0 || yi>M-1 ||
			si<0 || si>S-1 )
		{
			mexPrintf("Dropping %d: x=%d,y=%d,s=%d\n",p,xi,yi,si);
			continue;
		}

		// Each sample added to the histogram is weighted by its gradient magnitude
		// and by a Gaussian-weighted circular window with a \sigma that is 1.5 times
		// that of the scale of the keypoint.

		double sigma = sigma0*factor*pow(2,((double)s)/S0);
		int Gau_window = (int)floor(3.0*sigma);

		memset(Hist_pt,0,sizeof(Hist_pt));

		const double* pt = DoGs_ptr + yi*Y_OFFSET + xi*X_OFFSET + si*S_OFFSET;

#define at(dx,dy) (*(pt + (dx)*X_OFFSET + (dy)*Y_OFFSET))

		// in the Gaussian window
		for(int i=max(-Gau_window,1-xi);i<=min(Gau_window,N-2-xi);i++)
		{
			for (int j=max(-Gau_window,1-yi); j<=min(Gau_window,M-2-yi); j++)
			{
				double Dx = 0.5 * (at(i+1,j) - at(i-1,j));
				double Dy = 0.5 * (at(i,j+1) - at(i,j-1));
				double dx = ((double)(xi+i)) -x;
				double dy = ((double)(yi+j)) -y;

				if(dx*dx+dy*dy < Gau_window*Gau_window+0.5)
				{
					double theta = fmod(atan2(Dy,Dx)+2*PI,2*PI);
					int bin = (int) (BINS_SIZE*theta/(2*PI));
					double magnitude = sqrt(Dx*Dx+Dy*Dy);
					double Gau = exp(-(dx*dx+dy*dy)/(2*sigma*sigma));
					Hist_pt[bin] += magnitude*Gau; 
				}
			}
		}
		
		//////////////////////////////////////////////////////////////////////////

		// Smoothing histogram

		for (int i=0;i<6;i++)
		{
			double prev;
			prev = Hist_pt[BINS_SIZE-1];
			for(int j=0;j<BINS_SIZE;j++)
			{
				double val = (prev + Hist_pt[j]+Hist_pt[(j+1) % BINS_SIZE]) / 3.0;
				prev = Hist_pt[j];
				Hist_pt[j] = val;
			}
		}

		//////////////////////////////////////////////////////////////////////////

		// find highest peaks
		double max_hist = -1;
		for (int i=0;i<BINS_SIZE;i++)
			max_hist = max(max_hist,Hist_pt[i]);

		// any other local peak that is within 80% of the highest peak
		// is used to also create a keypoint with that orientation
		for (int i=0;i<BINS_SIZE;i++)
		{
			double curr  = Hist_pt[i];
			double left  = Hist_pt[(i-1+BINS_SIZE) % BINS_SIZE];
			double right = Hist_pt[(i+1+BINS_SIZE) % BINS_SIZE];

			if(curr > 0.8*max_hist && curr > left && curr > right)
			{
				double delta = 0.5*(right-left) / (2*curr-right-left);
				double val = 2*PI*(i+0.5+delta)/BINS_SIZE;

				// allocate memory dynamically
				if(buf_it == buf_tail)
				{
					int offset = buf_it - buf_head;
					buf_size += 4*max(1,BINS_SIZE/16);
					buf_head = (double*) mxRealloc(buf_head,buf_size*sizeof(double));
					buf_tail = buf_head + buf_size;
					buf_it = buf_head + offset;
				}
				// only about 15% of points are assigned multiple orientations
				*buf_it++ = x;
				*buf_it++ = y;
				*buf_it++ = s;
				*buf_it++ = val;
			}
		}

	}
	// Copy the result into an array.
	double* result;
	int NL = (buf_it - buf_head)/4 ;
	out[OUT_Q] = mxCreateDoubleMatrix(4, NL, mxREAL) ;
	result = mxGetPr(out[OUT_Q]);
	memcpy(result, buf_head, sizeof(double) * 4* NL) ;
	mxFree(buf_head) ;
}