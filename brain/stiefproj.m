function [x] = stiefproj(C)
[n,p] = size(C);
 
% Create the problem structure.
manifold = stiefelfactory(n, p);
problem.M = manifold;
 
% Define the problem cost function and its Euclidean gradient.
problem.cost  = @(x) trace((x-C)*(x-C)');
problem.egrad = @(x) 2*(x-C);
problem.ehess = @(x, xdot) 2*xdot;

x = trustregions(problem);

end