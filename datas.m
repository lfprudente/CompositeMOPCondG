function [dimA,A,b] = datas(n,m,delta)

dimA(1:m) = 2*n;

for ind = 1:m
    B = rand(n);
    A{ind} = [B; -B];
    b{ind}   = delta * ones(2*n,1);
end