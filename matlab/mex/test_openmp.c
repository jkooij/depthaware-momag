// Compile with 
//  >>> mex -v COMPFLAGS="$COMPFLAGS /openmp" test_openmp.c
// See http://nl.mathworks.com/matlabcentral/newsreader/view_thread/238799

#include "mex.h"
#include <omp.h>
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    omp_set_dynamic(1);
    // omp_set_num_threads(6);
    
    #pragma omp parallel
    {
        int i;
        mexPrintf("Max num threads %d.\n", omp_get_max_threads());

        #pragma omp for
        for ( i = 0; i < 10; i++)
        {
            mexPrintf("Num threads %d, thread ID %d.\n", omp_get_num_threads(), omp_get_thread_num());
        }
    }
}
