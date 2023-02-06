function [g,info] = evalg(n,x,ind,A,b)

info = 0;

options.Algorithm = 'dual-simplex';
options.Display = 'off';

[xopt,fmin,flag] = linprog(-x,A{ind},b{ind},[],[],-inf(n,1),inf(n,1),options);

if ( flag ~= 1 )
    info = -1;
    disp('ERROR!!')
    return
end

g = -fmin;