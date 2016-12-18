% mex -v hello.c CFLAGS="\$CFLAGS -fopenmp -lmwblas -lmwlapack" LDFLAGS="\$LDFLAGS -fopenmp -lmwblas -lmwlapack" -largeArrayDims

% mex -v hello.c -largeArrayDims -lmwblas -lmwlapack CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
% % mex -v sumofnorms.c -largeArrayDims -lmwblas -lmwlapack CFLAGS="\$CFLAGS -fopenmp -O3" LDFLAGS="\$LDFLAGS -fopenmp"
% 
% 
% setenv('LD_RUN_PATH', '/homes/mb2215/mkl')                   
% mex -v sumofnorms.c -largeArrayDims -lsingle_mkl_ilp64 ...
%     CFLAGS="\$CFLAGS -I /opt/intel/mkl/include -fopenmp" ...
%     LDFLAGS="\$LDFLAGS -L /homes/mb2215/mkl -fopenmp" ...
%     -DMKL_ILP64
% 
% mex -v maxerror.c -largeArrayDims -lmwblas -lmwlapack CFLAGS="\$CFLAGS -fopenmp -std=c99" LDFLAGS="\$LDFLAGS -fopenmp -std=c99"

mex -v quadsum.c -largeArrayDims CFLAGS="\$CFLAGS -fopenmp -std=c99" LDFLAGS="\$LDFLAGS -fopenmp -std=c99"