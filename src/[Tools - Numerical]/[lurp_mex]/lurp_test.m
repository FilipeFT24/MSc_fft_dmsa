% This function tests the routine lurp which implements Gaussian 
% elimination with rook pivoting using efficient code.
%
% To test standard rook pivoting use:
%
%    lurp_test
%
% To test threshold rook pivoting (type "help lurp" for information
% about the threshold), use:
%
%    lurp_test(tol)
% 
%  To test matrices of dimension m by m, m by 2*m and 2*m by m (the
%  default value is m = 200) use:
%
%    lurp_test(tol,m)
%
%  To return a vector of errors that should be of order of magnitude
%  machine epsilon ( typically smaller than m * eps ) use
%
%    errors=lurp_test, errors=lurp_test(tol), errors=lurp_test(tol,m)
%
% For large values of tol some of the tests may fail since
% threshold rook pivoting can introduce errors for large values of 
% tol.  Values of tol smaller than 10, including all negative values,
% should pass all tests.

%  created or modified by L. Foster, 7-10-2012, copr.
function errors = lurp_test(tol,m)

format compact
errors0 = [];
if nargin <=1
    m = 200;
end
m0 = m;
m = m0;
n = m0;
disp(' ')
disp(['m = ',num2str(m),', n = ',num2str(n)])
A = randn(m,n);
if nargin == 0
    e = test_matrix(A);
else
    e = test_matrix(A,tol);
end
errors0 = [errors0, e];

m = 2*m0;
n = m0;
disp(' ')
disp(['m = ',num2str(m),', n = ',num2str(n)])
A = randn(m,n);
if nargin == 0
    e = test_matrix(A);
else
    e = test_matrix(A,tol);
end
errors0 = [errors0, e];

m = m0;
n = 2*m0;
disp(' ')
disp(['m = ',num2str(m),', n = ',num2str(n)])
A = randn(m,n);
if nargin == 0
    e = test_matrix(A);
else
    e = test_matrix(A,tol);
end
errors0 = [errors0, e];

m = m0;
n = m0;
disp(' ')
if m > 5
   disp(['m = ',num2str(m),', n = ',num2str(n),' singular matrix '])
else
   disp(['m = ',num2str(m),', n = ',num2str(n)])
end
B=randn(m,5);
if m > 4
   B(3:5,:)=0;
end
C=randn(5,n);
if m > 4
   C(:,1:3)=0;
end
A=B*C;
if m > 4
   A(m-4:m,n-4:n)=A(m-4:m,n-4:n)+1.e-5*randn(5,5);
end
if nargin == 0
    e = test_matrix(A);
else
    e = test_matrix(A,tol);
end
errors0 = [errors0, e];

if nargout == 0
    disp(' ')
    if  ( max( errors0 ) < m0*eps )
        disp('lurp passed all tests')
    else
        nfail = sum( errors0 >= m0*eps );
        ntry = length(errors0);
        disp(['lurp failed ',int2str(nfail),' tests of ',int2str(ntry)])
    end
end

if nargout > 0
    errors = errors0;
end
 
end

function errors = test_matrix(A,tol)
errors = [];
% test with two outputs
if nargin == 1
   [L,U]=lurp(A);
else
   [L,U]=lurp(A,tol);
end
normA = norm(A);
e = norm(A-L*U)/normA;
disp(['with two outputs: error = ',num2str(e)])
errors = [errors, e];
% test with three outputs
if nargin == 1
   [L,U,P]=lurp(A);
else
   [L,U,P]=lurp(A,tol);
end
e = norm(P*A-L*U)/normA;
errors = [errors, e];
disp(['with three outputs: error = ',num2str(e)])
% test with vector for permutation
if nargin == 1
   [L,U,p]=lurp(A,'vector');
else
   [L,U,p]=lurp(A,'vector',tol);
end
e = norm(A(p,:)-L*U)/normA;
errors = [errors, e];
disp(['with permutation stored as a vector: error = ',num2str(e)])
% test with four outputs
if nargin == 1
   [L,U,P,Q]=lurp(A);
else
   [L,U,P,Q]=lurp(A,tol);
end
e = norm(P*A*Q-L*U)/normA;
errors = [errors, e];
disp(['with four outputs: error = ',num2str(e)])
% test with vectors for permutations
if nargin == 1
   [L,U,p,q]=lurp(A,'vector');
else
   [L,U,p,q]=lurp(A,'vector',tol);
end
e = norm(A(p,q)-L*U)/normA;
errors = [errors, e]; 
disp(['with permutations stored as vectors: error = ',num2str(e)])
end

