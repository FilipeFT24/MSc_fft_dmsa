% The LURP package has efficient Fortran code to compute
% LU factorizations of dense, unsymmetric matrices using
% Gaussian elimination with rook pivoting. Also the package 
% includes a MATLAB mex interface so that the routines
% can be called directly from MATLAB as well as MATLAB
% code to install and test the package.
%
% The LURP package was created by L. Foster, Department of
% Mathematics, San Jose State University. 7-31-2012, copr.
% Use under the GNU General Public License (GPL) guidelines.
%
% FILES:
%
% lurp_install.m     : installs the LURP package
% lurp.m             : "help lurp" describes the use of LURP
% lurp.F             : the Mex gateway routine for LURP package
% lurpf.F            : a Mex "computational routine" called by lurp.F
% dgerpf.F           : Fortran code implementing rook pivoting
% dgerp2.F           : unblocked rook pivoting, called by DGERPF
% dgerp3.F           : blocked, partial factorization, called by DGERPF
% dlaswq.F           : permutes columns of a matrix, called by DGERPF
% ilaenv.F           : Fortran version of LAPACK code that returns
%                    : the block size. From LAPACK 3.4. Included since
%                    : the Fortran version returns a good block size
%                    : whereas the version of ILAENV that comes with
%                    : MATLAB's LAPACK library may not.
% lurp_test.m        : quick test of the compiled LURP
% lurp_test_SJ.m     : test of LURP using matrices from  the San Jose
%                    : State University Sparse Matrix Database
% lurp_test_UF.m     : test of LURP using matrices from the University
%                    : of Florida Sparse Matrix Collection
% readme.txt         : describes the installation and use of the package
% Contents.m         : this file


