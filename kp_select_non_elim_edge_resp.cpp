
#include"mex.h"

#include<stdlib.h>
#include<string.h>

//////////////////////////////////////////
#include<math.h>

const int MAX_ITER = 5 ;

#define abs(a)       (((a)>0)?(a):(-a))

void mexFunction(int nout, mxArray *out[], int nin, const mxArray *in[])
{      
  const double* Loc_Max_ptr ;
  const double* DoGs_ptr ;
  const int* dims;
  double threshold = 0.01 ; //0.01 the extrema not less than threshold.
  double r = 10.0 ;  // Eliminating edge responses, threshold.
  double* result ;
  enum {IN_LOC=0,IN_DOGS,IN_THRESHOLD,IN_R} ;
  enum {OUT_Q=0} ;  
  if(nin >= 3) 
  {
    threshold = *mxGetPr(in[IN_THRESHOLD]) ;
  }
  if(nin >= 4) 
  {
    r = *mxGetPr(in[IN_R]);
  }  
  // DoGs 
  int M,N,S;
  dims = mxGetDimensions(in[IN_DOGS]);
  M = dims[0] ; //y range [1,M]
  N = dims[1] ; //x range [1,N]
  S = dims[2] ; //scale   [1,S]
 if(S < 3 || M < 3 || N < 3) {
    mexErrMsgTxt("All dimensions of DOG must be not less than 3.") ;
  }
  DoGs_ptr = mxGetPr(in[IN_DOGS]) ; // head pointer of DoGs vectors (M*N*S)
  
  // Local Max Points
  int LOCAL_MAX_NUM = mxGetN(in[IN_LOC]) ; //the number of Local Maximum points .
  Loc_Max_ptr = mxGetPr(in[IN_LOC]) ; // head pointer of Local Max vector (3*LOCAL_MAX_NUM)

  // If the input array is empty, then output an empty array as well.
  if( LOCAL_MAX_NUM == 0) {
    out[OUT_Q] = mxDuplicateArray(in[IN_LOC]) ;
    return ;
  }
 
  ///////////////////////////////////////////////////////////////////////

    double* buffer = (double*) mxMalloc(LOCAL_MAX_NUM*3*sizeof(double)) ;
    double* buffer_iterator = buffer ;
	// (M,N,S)
    const int Y_OFFSET = 1;
    const int X_OFFSET = M;
    const int S_OFFSET = M*N ;

    for(int p = 0 ; p < LOCAL_MAX_NUM ; ++p) {

      int x = ((int)*Loc_Max_ptr++) ;
      int y = ((int)*Loc_Max_ptr++) ;
      int s = ((int)*Loc_Max_ptr++) ;
      double _offset[3] ;

#define at(dx,dy,ds) (*(pt + (dx)*X_OFFSET + (dy)*Y_OFFSET + (ds)*S_OFFSET))

	  const double* pt = DoGs_ptr + y*Y_OFFSET + x*X_OFFSET + s*S_OFFSET;
	  double Dx=0,Dy=0,Ds=0,Dxx=0,Dyy=0,Dss=0,Dxy=0,Dxs=0,Dys=0 ;
      int dx = 0;
      int dy = 0;
	  // find the sample point that th offset (dx,dy,ds)
	  for(int iter = 0 ; iter < MAX_ITER ; ++iter) 
	  {
		  double H[3*3] ;
		  x += dx ;
		  y += dy ;
		  pt = DoGs_ptr + y*Y_OFFSET + x*X_OFFSET + s*S_OFFSET ;		  

		  // Compute the gradient.
		  Dx = 0.5 * (at(+1,0,0) - at(-1,0,0)) ;
		  Dy = 0.5 * (at(0,+1,0) - at(0,-1,0));
		  Ds = 0.5 * (at(0,0,+1) - at(0,0,-1)) ;

		  // Compute the Hessian.
		  Dxx = (at(+1,0,0) + at(-1,0,0) - 2.0 * at(0,0,0)) ;
		  Dyy = (at(0,+1,0) + at(0,-1,0) - 2.0 * at(0,0,0)) ;
		  Dss = (at(0,0,+1) + at(0,0,-1) - 2.0 * at(0,0,0)) ;

		  Dxy = 0.25 * ( at(+1,+1,0) + at(-1,-1,0) - at(-1,+1,0) - at(+1,-1,0) ) ;
		  Dxs = 0.25 * ( at(+1,0,+1) + at(-1,0,-1) - at(-1,0,+1) - at(+1,0,-1) ) ;
		  Dys = 0.25 * ( at(0,+1,+1) + at(0,-1,-1) - at(0,-1,+1) - at(0,+1,-1) ) ;

		  //Hessian matrix
		  H[0+3*0] = Dxx;
		  H[1+3*1] = Dyy;
		  H[2+3*2] = Dss;
		  H[0+3*1] = H[1+3*0] = Dxy;
		  H[0+3*2] = H[2+3*0] = Dxs;
		  H[1+3*2] = H[2+3*1] = Dys;

		  // solve linear system
		  _offset[0] = - Dx ;
		  _offset[1] = - Dy ;
		  _offset[2] = - Ds ;

		  // Gauss elimination //row reduction
		  int i,j,ii,jj;
		  for(j = 0 ; j < 3 ; ++j) {
			  double Max_h		= 0 ;
			  double Max_abs_h	= 0 ;
			  int	 max_idx	= -1 ;			  

			  // look for the maximally stable pivot
			  for (i = j ; i < 3 ; ++i) {			  
				  if ( abs(H[i+3*j]) > Max_abs_h) {
					  Max_h = H[i+3*j];
					  Max_abs_h = abs(H[i+3*j]) ;
					  max_idx = i;
				  }
			  }
			  // if singular give up
			  if (Max_abs_h < 1e-10f) {
				  _offset[0] = 0 ;
				  _offset[1] = 0 ;
				  _offset[2] = 0 ;
				  break ;
			  }
			  double tmp;
			  // swap j-th row with i-th row and normalize j-th row
			  for(jj = j ; jj < 3 ; ++jj) {
				  tmp = H[max_idx+3*jj] ; H[max_idx+3*jj] = H[j+3*jj] ; H[j+3*jj] = tmp ;
				  H[j+3*jj] /= Max_h;
			  }
			  tmp = _offset[j] ; _offset[j] = _offset[max_idx] ; _offset[max_idx] = tmp ;
			  _offset[j] /= Max_h ;

			  // elimination
			  for (ii = j+1 ; ii < 3 ; ++ii) {
				  double x = H[ii+3*j] ;
				  for (jj = j ; jj < 3 ; ++jj) {
					  H[ii+3*jj] -= x * H[j+3*jj] ;
				  }
				  _offset[ii] -= x * _offset[j] ;
			  }
		  }

		  // backward substitution
		  for (i = 2 ; i > 0 ; --i) {
			  double x = _offset[i] ;
			  for (ii = i-1 ; ii >= 0 ; --ii) {
				  _offset[ii] -= x * H[ii+3*i] ;
			  }
		  }

		  // If the translation of the keypoint is big, move the keypoint
		  // and re-iterate the computation. Otherwise we are all set.
		  //
		  // 0.6 -> 0.5
		  double thes = 0.6;
		  dx= ((_offset[0] >  thes && x < N-2) ?  1 : 0 )
			+ ((_offset[0] < -thes && x > 1  ) ? -1 : 0 ) ;

		  dy= ((_offset[1] >  thes && y < M-2) ?  1 : 0 )
			+ ((_offset[1] < -thes && y > 1  ) ? -1 : 0 ) ;

		  if( dx == 0 && dy == 0 ) break ;

	  }

	  // the extrema value
	  double val = at(0,0,0) + 0.5 * (Dx * _offset[0] + Dy * _offset[1] + Ds * _offset[2]) ;

	  double score = (Dxx+Dyy)*(Dxx+Dyy) / (Dxx*Dyy - Dxy*Dxy) ;
	  double xn = x + _offset[0] ; 
	  double yn = y + _offset[1] ;
	  double sn = s + _offset[2] ;

	  if( fabs(val) > threshold &&
		  fabs(_offset[0]) < 1.5 && 
		  fabs(_offset[1]) < 1.5 && 
		  fabs(_offset[2]) < 1.5 &&
		  xn >= 0 && xn <= N-1 && 
		  yn >= 0 && yn <= M-1 && 
		  sn >= 0  && sn <= S-1) 
	  {
		  *buffer_iterator++ = xn;
		  *buffer_iterator++ = yn;
		  *buffer_iterator++ = sn;
	  }
	}

	// Copy the result into an array.
	{
		int NL = (buffer_iterator - buffer)/3 ;
		out[OUT_Q] = mxCreateDoubleMatrix(3, NL, mxREAL) ;
		result = mxGetPr(out[OUT_Q]);
		memcpy(result, buffer, sizeof(double) * 3 * NL) ;
	}
	mxFree(buffer) ;
}


