function [theta,p,flag] = subproblemCondG(n,m,l,u,x,JH,G,dimA,A,b)

flag = -1;

sdimA = sum(dimA);

c = [1; zeros(n+sdimA,1)];

index(1) = n + 1;
for ind = 2:m+1
    index(ind) = index(ind-1) + dimA(ind-1);
end

Asubprob   = zeros(m,1+n+sdimA);
Aeqsubprob = zeros(n*m,1+n+sdimA);

for ind = 1:m
    Asubprob(ind,1:1+n) = [-1 JH(ind,:)];
    Asubprob(ind,index(ind)+1:index(ind+1)) = b{ind};
    Aeqsubprob((ind-1)*n+1:ind*n,2:n+1) = -eye(n);
    Aeqsubprob((ind-1)*n+1:ind*n,index(ind)+1:index(ind+1)) = A{ind}';
end

for ind = 1:m
    bsubprob(ind) = G(ind) + dot(JH(ind,:),x); 
end
bsubprob = bsubprob';

beqsubprob(1:n*m) = 0;
beqsubprob = beqsubprob';

lb = [-inf; l; zeros(sdimA,1)];
ub = [ inf; u; inf(sdimA,1)  ];

options.Algorithm = 'dual-simplex';
options.Display = 'off';

[xopt, fmin, info] = linprog(c, Asubprob, bsubprob, Aeqsubprob, beqsubprob,lb, ub,options);

if ( info == 1 )
    flag = 0;
end

theta = xopt(1);
p = xopt(2:n+1);