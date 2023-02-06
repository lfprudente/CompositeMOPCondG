function [x,info,iter,nfev,ngev,time] = CondG(n,m,l,u,x,dimA,A,b)

tic;

tol  = 10^(-4);
ftol = 10^(-4);

itermax = 200;

% Counters

nfev = 0;
ngev = 0;
iter = 0;

% Evaluate F

for i = 1:m
    H(i) = evalh(n,x,i);
    G(i) = evalg(n,x,i,A,b);
end
F = H + G;
nfev = nfev + 1;

fprintf('----------------------------------------------------------------------\n')
fprintf('   Proximal Gradient Method for Composite Multiojective Optimization  \n')
fprintf('----------------------------------------------------------------------\n')
fprintf('Number of variables  : %i \n',n)
fprintf('Number of objectives : %i \n',m)
fprintf('Optimality tolerance: %.0e \n\n',tol)

% -----------
% Main Loop
% -----------

while(1)

    % Compute the Jacobian of H

    for i = 1:m
        JH(i,:) = evalgradh(n,x,i);
    end
    ngev = ngev + 1;

    % Solve the subproblem

    [theta,p,flagIS] = subproblemProxGrad(n,m,l,u,x,JH,G,dimA,A,b);

    % Print information

    if ( mod(iter,10) == 0 )
        fprintf('\n')
        fprintf('%-5s   %-8s   %-8s%-8s \n','it','|theta|','IS','LS')
    end
    if ( iter == 0 )
        fprintf('%5d   %8.2e   %-8i %-8s\n',iter,abs(theta),flagIS,'-')
    else
        fprintf('%5d   %8.2e   %-8i %-8i\n',iter,abs(theta),flagIS,flagLS)
    end
    
    % -----------
    % Stopping criteria
    % -----------

    % Test optimality 
    
     if ( iter > 0 && norm( x - xprev, Inf ) / max(1 , norm( xprev, Inf ) ) <= tol )

        if ( abs(theta) <= tol )
            info = 0;

            % Stop timing 

            time = toc;

            % Print information

            fprintf('\n')
            fprintf('Solution was found.\n')
            fprintf('Number of function evaluations: %i\n',nfev)
            fprintf('Number of gradient evaluations: %i\n',ngev)
            fprintf('CPU time(s)                   : %.1f \n',time)

            return
        end
    
     end
        
    % Test the maximum number of iterations

    if ( iter >= itermax )
        info = 1;

        % Stop timing 

        time = toc;

        % Print information

        fprintf('\n')
        fprintf('The number of maximum iterations was reached.\n')
        fprintf('CPU time(s): %.1f \n',time)

        return
    end

    % -----------
    % Iterate
    % -----------

    iter = iter + 1;

    % Define the search direction

    d = p - x;

    % Compute the step size

    [stp,H,G,nfevLS,flagLS] = armijo(n,m,x,d,F,JH,theta,ftol,A,b);
    nfev = nfev + nfevLS;

    %stp = min(1,abs(theta)/(2* norm(p-x)^2));

    % Update x
    
     xprev = x;

    x = x + stp * d;

    % Evaluate F
 
    F = H + G;

end

% -----------
% End of Main Loop
% -----------
