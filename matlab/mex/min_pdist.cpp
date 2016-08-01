#include <mex.h>
#include <cmath>
#include <limits>

// #include "mex_eigen_utils.hpp"
#include "StringBuilder.h"

void mexPrintString(const std::string &s)
{
    mexPrintf(s.c_str());
}

#define MEX_PRINT_STRING(x) mexPrintString(UTIL_STR(x));
#define DEBUG_SHOW(x) MEX_PRINT_STRING(#x << " = " << (x) << "\n");

// X : nx x d
// Y : ny x d
// [mindists, minidxs] = min_pdist_thresh(X, Y, dthresh)
//   minidxs: ny x 1
//   mindists: ny x 1

// compute squared euclidean distance
double eucldist2(double *px, double *py, int dims, int deltapx, int deltapy)
{
    int d;
    double dist = 0.;
    
    for(d = 0; d < dims; ++d) {
        dist += (*px - *py) * (*px - *py);
        px += deltapx;
        py += deltapy;
    }
    
    return dist;
}

double eucldist2(double *px, double *py, int dims, int deltapx, int deltapy, double bestdist)
{
    int d;
    double dist = 0.;
    
    for(d = 0; (d < dims) && (dist < bestdist); ++d) {
        dist += (*px - *py) * (*px - *py);
        px += deltapx;
        py += deltapy;
    }
    
    return dist;
}

double cityblockdist2(double *px, double *py, int dims, int deltapx, int deltapy)
{
    int d;
    double dist = 0.;
    
    for(d = 0; d < dims; ++d) {
        dist += std::abs(*px - *py);
        px += deltapx;
        py += deltapy;
    }
    
    return dist * dist;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // declare variables
    const mxArray *mX, *mY;
    mxArray *mMinDists, *mMinIdxs;
    
    const mwSize *dimX, *dimY;
    int nx, ny, D;
    double *pX, *pY, *pMinDists, *pMinIdxs;
    
    double dthresh = 180. * 180.; // distance threshold

    // check for proper number of arguments
    if (nrhs != 2) { 
        mexErrMsgTxt("Two input argument required."); 
    }  
    if (nlhs != 2) {
        mexErrMsgTxt("Two outputs required."); 
    } 

    // setup inputs
    mX = prhs[0]; mY = prhs[1];
    dimX = mxGetDimensions(mX);
    dimY = mxGetDimensions(mY);
    
    nx = (int)dimX[0];
    ny = (int)dimY[0];
    D = (int)dimX[1];
    
    if (D != dimY[1]) {
        mexErrMsgTxt("Different number of arguments"); 
    }

    pX = mxGetPr(mX);
    pY = mxGetPr(mY);
    
    // setup outputs
    mMinIdxs = plhs[0] = mxCreateDoubleMatrix(ny, 1, mxREAL);
    mMinDists = plhs[1] = mxCreateDoubleMatrix(ny, 1, mxREAL);
    
    pMinDists = mxGetPr(mMinDists);
    pMinIdxs = mxGetPr(mMinIdxs);
    
    // test all distances
    int offset = 0;
    //evals = 0;
    for(int yj = 0; yj < ny; yj++)
    {
        int bestidx = std::numeric_limits<double>::max();
        double bestdist = std::numeric_limits<double>::max();
        //evals = 0;
        for(int xj2 = 0; xj2 < nx; xj2++)
        {
            int xj = (xj2 + offset) % nx;
            
            
            if (xj > 0) {
                double dist = cityblockdist2(pX+xj, pY+yj, D, nx, ny);
                if (dist > dthresh) {
                    continue;
                }
            }
            
            
            double dist = eucldist2(pX+xj, pY+yj, D, nx, ny, bestdist);
            //MEX_PRINT_STRING("xj = " << xj << "\tbestdist = " << bestdist << "\tdist = " << dist << "\n");
            
            if (dist < bestdist)
            {
                bestdist = dist;
                bestidx = xj;
            }            
        }
        
        offset = bestidx;
        //DEBUG_SHOW(bestidx);
        //DEBUG_SHOW(evals)
        
        
        *pMinDists = std::sqrt(bestdist);
        *pMinIdxs = bestidx+1;
        pMinDists++;
        pMinIdxs++;
    }
    
    //DEBUG_SHOW(evals);
}
