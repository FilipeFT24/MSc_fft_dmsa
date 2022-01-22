addpath(genpath('[Tools - Numerical]'));

  fd = @(p) p(:,1).^2/2^2 + p(:,2).^2/1^2 - 1;
  fh = @(p) ones(size(p,1),1);
  [p,t] = distmesh2d( fd, fh, 0.2, [-2,-1;2,1],[]);