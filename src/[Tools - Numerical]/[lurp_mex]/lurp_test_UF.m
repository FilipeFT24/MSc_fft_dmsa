% This function tests the routine lurp which implements Gaussian 
% elimination with rook pivoting using efficient code.  The test matrices
% come from the University of Florida sparse matrix collection. The
% matrices are converted to dense matrices using the Matlab full command.
%
% To test standard rook pivoting on the 200 smallest dimension matrices
% from the database:
%
%    lurp_test_UF
%
% This function requires the installation of UFget which is available at
% http://www.cise.ufl.edu/research/sparse/mat/UFget.html .
%
% For each matrix lurp is tested with two outputs, three outputs with
% permutations stored as matrices, three outputs with permutations
% stored a vectors, and four outputs with permutations stored as
% matrices and as vectors.
%
% To test treshold rook pivoting use:
%
%    lurp_test_UF(tol)
% 
% To test the nsamples smallest dimension matrices from the database use:
%
%    lurp_test_UF(tol,nsamples)
%
% For large values of tol some of the tests may fail since
% threshold rook pivoting can introduce errors for large values of 
% tol.  Values of tol smaller than 10, including all negative values,
% should pass all tests.

% If lurp_test_UF has an output parameter the parameter contains errors
% that should be of order of magnitude machine epsilon (typically smaller
% than max(size(A)) * eps ). With two output parameters the second
% parameter is a vector of the indices in the UF database corresponding
% to each component of the percent error.
%
% Only real parts of complex matrices are used for the tests.

%  created or modified by L. Foster, 7-10-2012, copr.

function [errors, UFindices] = lurp_test_UF(tol,nsamples)

location = which('UFget');
if isempty(location)
    disp('UFget is required and is not installed.')
    disp('To install the UFget tools see')
    disp('http://www.cise.ufl.edu/research/sparse/mat/UFget.html')
    return
end
    
if nargin <= 1
   nsamples=200;
end
format compact
format short e
errors0 = [];

index = UFget;
dim = max(index.nrows,index.ncols);
[dim_sort,indexs]=sort(dim);
indexs = indexs(1:nsamples);
dim_sort=dim_sort(1:nsamples);
cnt = 0;
dims = [];
ifail = [];
disp(' ')
for i = indexs
    cnt = cnt + 1;
    Problem = UFget(i,index);
    A=Problem.A;
    A = real(A);   % no complex samples
    [m,n]=size(A);
    disp(['Matrix ',int2str(cnt),' of ',int2str(nsamples),':', ...
        ' size(A) = ',int2str(m),' ',int2str(n)])
    A = full(A);
    if nargin == 0
        e = test_matrix_inf(A);
    else
        e = test_matrix_inf(A,tol);
    end
    errors0 = [errors0, e];
    dims_new = dim_sort(cnt)*ones(size(e));
    if ( sum ( e >= dims_new*eps ) ~= 0)
        ifail = [ ifail, i ];
    end
    dims = [ dims, dims_new ];

end
if nargout == 0
    disp(' ')
    if  ( sum( errors0  >= dims*eps ) == 0)
        disp('lurp passed all tests')
    else
        nfail = sum( errors0 >= dims*eps );
        ntry = length(errors0);
        disp(['lurp failed ',int2str(nfail),' tests of ',int2str(ntry)])
        disp('at the matrices with indices in the UF database = ')
        disp( int2str(ifail) )
    end
end
if nargout > 0
    errors = errors0;
    UFindices = [indexs; indexs; indexs; indexs; indexs];
    UFindices = UFindices(:)';
end

end

function errors = test_matrix_inf(A,tol)
errors = [];
% test with two outputs
if nargin == 1
   [L,U]=lurp(A);
else
   [L,U]=lurp(A,tol);
end
normA = norm(A,inf);
e = norm(A-L*U,inf)/normA;
errors = [errors, e];
% test with three outputs
if nargin == 1
   [L,U,P]=lurp(A);
else
   [L,U,P]=lurp(A,tol);
end
e = norm(P*A-L*U,inf)/normA;
errors = [errors, e];
if nargin == 1
   [L,U,p]=lurp(A,'vector');
else
   [L,U,p]=lurp(A,'vector',tol);
end
e = norm(A(p,:)-L*U,inf)/normA;
errors = [errors, e];
if nargin == 1
   [L,U,P,Q]=lurp(A);
else
   [L,U,P,Q]=lurp(A,tol);
end
e = norm(P*A*Q-L*U,inf)/normA;
errors = [errors, e];
if nargin == 1
   [L,U,p,q]=lurp(A,'vector');
else
   [L,U,p,q]=lurp(A,'vector',tol);
end
e = norm(A(p,q)-L*U,inf)/normA;
errors = [errors, e]; 
end

