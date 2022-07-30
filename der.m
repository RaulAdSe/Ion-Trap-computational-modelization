function d = der(fun,x)
h = 10^-5;
d = (fun(x+h)-fun(x))/h;
end
