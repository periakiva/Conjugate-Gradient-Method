function [u, niter, residue] = solveCG(A, f, tol,maxiter)
%
%    Input parameters: 
%           A : Symmetric, positive definite NxN matrix 
%           f : Right-hand side Nx1 column vector 
%         tol : relative residual error tolerance for break
%               condition 
%     maxiter : Maximum number of iterations to perform
%
%    Output parameters:
%           u : Nx1 solution vector
%       niter : Number of iterations performed
%        flag : 1 if convergence criteria specified by TOL could
%               not be fulfilled within the specified maximum
%               number of iterations, 0 otherwise (= iteration
%               successful).

s=f;
u = s;         % Set u_0 to the start vector s
r = f - A*s;   % Compute first residuum
p = r;         
rho = r'*r;
niter = 0;     % Init counter for number of iterations
flag = 0;      % Init break flag
residue=[];

% Compute norm of right-hand side to take relative residuum as
% break condition.
normf = norm(f);
if normf < eps  % if the norm is very close to zero, take the
                % absolute residuum instead as break condition
                % ( norm(r) > tol ), since the relative
                % residuum will not work (division by zero).
  warning(['norm(f) is very close to zero, taking absolute residuum' ... 
					 ' as break condition.']);
	normf = 1;
end

while (norm(r)/normf > tol)   % Test break condition
	a = A*p;
    
	alpha = rho/(a'*p);
	u = u + alpha*p;
	r = r - alpha*a;
    res = mean(r);
	rho_new = r'*r;
	p = r + rho_new/rho * p;
	rho = rho_new;
    %pause(0.5)
    residue = [residue;rho];
	niter = niter + 1;
    %residue=[residue;norm(r)];
	if (niter == maxiter)         % if max. number of iterations
		flag = 1;                   % is reached, break.
		break
	end
end