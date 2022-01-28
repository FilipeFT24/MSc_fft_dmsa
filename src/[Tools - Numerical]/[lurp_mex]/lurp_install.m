function errors = lurp_install
% This function installs lurp, a package that implements Gaussian
% elimination with rook pivoting using code that minimizes, where possible,
% memory access for increased efficiency.
%
% To install the package:
%
%    download and uncompress the zip file containing the package
%    start MATLAB and move from inside MATLAB to the folder
%        containing the uncompressed files
%    set the proper compiler for MATLAB's mex utility by typing
%        mex -setup 
%        and following the instructions
%    type
%        lurp_install
%
% The installation is designed to work for Linux systems using the
% gfortran compiler or Windows systems using Intel's Fortran compiler.
% It is assumed that MATLAB has been configured so that one of these
% compilers is the default compiler for MATLAB's mex facility.  The
% compiler should be supported by MATLAB.  See
%    http://www.mathworks.com/support/compilers/R2012a/win64.html
%      or
%    http://www.mathworks.com/support/sysreq/previous_releases.html
% 
% The installation has been tested on a 64 bit Linux system running
% Scientific Linux, a variant of Redhat Linux, with MATLAB 7.12 and
% GNU gfortran 4.3.  It has also been tested with 64 bit and 32 bit
% MATLAB 7.9 on a system running Windows XP and Intel Fortran 10.1.
% The installation script may work in other environments but has not
% been tested.

% To install the package and return a vector of errors from tests of
% the installation type
%    errors = lurp_install
% The entries in errors should be of order of magnitude machine epsilon
% (smaller than 200*eps for the tests in lurp_install).

% Written by Leslie Foster, 7-23-2012, copr.

version_string = version;
% assume at least two decimal points in version number
idot = strfind(version_string,'.');
version_num = str2double(version_string(1:idot(2)-1));
version_int = str2double(version_string(1:idot(1)-1));
version_frac = str2double(version_string(idot(1)+1:idot(2)-1));
if ( ( (version_int > 7) || (version_int == 7) && (version_frac >= 8) ) ...
        || (version_num == 7.6) && strcmp(mexext,'mexw64') || ...
        (version_num == 7.6) && strcmp(mexext,'mexa64') )
   disp(['Compiling lurp in Matlab version ',version_string])
else
   disp('These mex utilities are only designed for Matlab 7.8 (32 or')
   disp('   64 bit) or later or 64 bit Matlab 7.6.  Modifications of the')
   disp('   mex commands may work with other versions of Matlab.')
   return
end

% for 32 bit Matlab 7.8 or later for Windows
if ( ( (version_int > 7) || (version_int == 7) && (version_frac >= 8) ) ...
        && strcmp(mexext,'mexw32') )
   blaslib = fullfile(matlabroot,'extern', 'lib', 'win32',...
       'microsoft', 'libmwblas.lib');
   lapacklib = fullfile(matlabroot,'extern', 'lib', 'win32', ...
       'microsoft', 'libmwlapack.lib');
   mex('-v','lurp.F','lurpf.F','dgerpf.F','dgerp2.F','dgerp3.F',...
       'dlaswq.F','ilaenv.F',lapacklib,blaslib);
end 
% for 32 bit Matlab 7.8 or later for Linux
if ( ( (version_int > 7) || (version_int == 7) && (version_frac >= 8) ) ...
        && strcmp(mexext,'mexa32') )
   mex -v lurp.F lurpf.F dgerpf.F dgerp2.F dgerp3.F dlaswq.F...
       ilaenv.F -lmwlapack -lmwblas
end 

% In Matlab 7.6 or 7.7 it appears best not to use largeArrayDims 
% when using the BLAS or LAPACK libraries since these libraries assume that
% the array indices are not for large arrays.
% For 64 bit Matlab 7.6 (or potentially 7.7) for Windows:
if ( ( (version_int == 7) && (version_frac < 8) ) ...
        && strcmp(mexext,'mexw64') )
   blaslib = fullfile(matlabroot,'extern', 'lib', 'win64', ...
       'microsoft', 'libmwblas.lib');
   lapacklib = fullfile(matlabroot,'extern', 'lib', 'win64', ...
       'microsoft', 'libmwlapack.lib');
   mex('-v','lurp.F','lurpf.F','dgerpf.F','dgerp2.F','dgerp3.F',...
       'dlaswq.F','ilaenv.F',lapacklib,blaslib);
end 
% For 64 bit Matlab 7.6 (or potentially 7.7) for Linux:
if ( ( (version_int == 7) && (version_frac < 8) ) ...
        && strcmp(mexext,'mexa64') )
   mex -v lurp.F lurpf.F dgerpf.F dgerp2.F dgerp3.F dlaswq.F ...
       ilaenv.F -lmwlapack -lmwblas
end 

% If one is accessing the BLAS or LAPACK library that comes with Matlab 7.8
% or later it is required that integer parameters in calls to LAPACK or
% BLAS routines be 64bit integers. The additional compiler
% flags below cause the ifort compiler in Windows and the gfortran
% compiler in Linux to allocate 64 bit integers for all integer
% variables and integer constants. Therefore to implement largeArrayDims 
% and to access the BLAS or LAPACK in Matlab 7.8 or later it appears that 
% one does not need to manually make the code changes  suggested at 
% http://www.mathworks.com/support/solutions/data/1-5C27B9.html?solution=1-5C27B9
% With these compilers switches it appears that the Fortran code, outside 
% the MEX gateway routine lurp.F, can be standard Fortran code.
% Matlab 7.8 or later with largeArrayDims:
if ( ( (version_int > 7) || (version_int == 7) && (version_frac >= 8) ) ...
        && strcmp(mexext,'mexw64') )
   blaslib = fullfile(matlabroot,'extern', 'lib', 'win64', ...
       'microsoft', 'libmwblas.lib');
   lapacklib = fullfile(matlabroot,'extern', 'lib', 'win64', ...
       'microsoft', 'libmwlapack.lib');
   mex('-v','-largeArrayDims',...
       'COMPFLAGS="$COMPFLAGS','/4I8','/intconstant"','lurp.F',...
       'lurpf.F','dgerpf.F','dgerp2.F','dgerp3.F','dlaswq.F',...
       'ilaenv.F',lapacklib,blaslib);
end
if ( ( (version_int > 7) || (version_int == 7) && (version_frac >= 8) ) ...
        && strcmp(mexext,'mexa64') )
   mex -v FFLAGS="\$FFLAGS -fdefault-integer-8 -O3" -largeArrayDims ...
       lurp.F lurpf.F dgerpf.F dgerp2.F dgerp3.F dlaswq.F ilaenv.F ...
       -lmwlapack -lmwblas
end

disp('Finished compiling.')
disp(' ')
disp('Press enter to test the lurp code.')
disp(' ')
pause
disp('STANDARD ROOK PIVOTING:')
errors1 = lurp_test;
disp(' ')
disp('THRESHOLD ROOK PIVOTING, TOL = 2:')
errors2 = lurp_test(2);
disp(' ')
disp('THRESHOLD ROOK PIVOTING, TOL = -2:')
errorsm2 = lurp_test(-2);
m0 = 200;  % the default value of min(size(A)) in lurp_test
errors0 = [errors1,errors2,errorsm2];
disp(' ')
if  ( max( errors0 ) < m0*eps )
    disp('lurp passed all tests')
else
    nfail = sum( errors0 >= m0*eps );
    ntry = length(errors0);
    disp(['lurp failed ',int2str(nfail),' tests of ',int2str(ntry)])
end

if nargout > 0
    errors = errors0;
end

file1 = mfilename('fullpath');
[path1,name,extension]= fileparts(file1);
addpath(path1);
disp(' ')
disp(['Path ',path1,' added to the MATLAB path.'])
disp(' ')
disp('If you wish later access to the lurp routine type savepath.')
disp(' ');
