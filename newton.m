function [xk,it] = newton (a,tol,itmax,fun)
% a: initial guess
% tol: tolerance
% itmax: maximum number of iterations
% fun: function to enter
xk=[a] ; it=0 ; dx=1e-4; tolk=10;
while it < itmax && tolk > tol
    fk=fun(xk(end)) ; fk2 = fun(xk(end)+dx);
    dfk = (fk2-fk)/dx;
    xk = [xk  xk(end) - fk/dfk] ;   
    tolk = abs(xk(end)-xk(end-1));
    it=it+1;
end



