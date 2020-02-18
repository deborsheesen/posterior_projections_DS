function [x, xnorm] = stiefproj(C, ind)
[n,p] = size(C);

if ind == 1
    xnorm = norm(C,'fro');
    C = C*(1/xnorm);
else
    xnorm = norm(C);
    C = C*(1/xnorm);
end
 
% Create the problem structure.
manifold = stiefelfactory(n, p);
problem.M = manifold;
 
% Define the problem cost function and its Euclidean gradient.
problem.cost  = @(x) trace((x-C)*(x-C)');
problem.egrad = @(x) 2*(x-C);
problem.ehess = @(x, xdot) 2*xdot;

x = trustregions(problem);

end