%function [mmsh, E_cubic] = take8(dim, indraw, Vi)
prova = zeros(2,2,2);
prova(:,:,1) = [1,2;3,4];
prova(:,:,2) = [5,6;7,8];
indraw = [1.5 1.5 1.5];
dim = [0,1]

mi = zeros(3,2);
for k = 1:3
    mi(k, 1) = ceil(indraw(k));
    mi(k, 2) = floor(indraw(k));
end

[ix, iy, iz] = meshgrid([mi(1,1) mi(1,2)], [mi(2,1) mi(2,2)], [mi(3,1) mi(3,2)]);
[mx, my, mz] = meshgrid([dim(mi(1,1)) dim(mi(1,2))], [dim(mi(2,1)) dim(mi(2,2))], [dim(mi(3,1)) dim(mi(3,2))]);
mmsh = {mx, my, mz};

H = dim(2) - dim(1);
V8 = zeros(2,2,2);  V_cubic = V8;
for k = 1:8
    V8(k) = Vi(ix(k), iy(k), iz(k));
end
V8 = permute(V8, [2 1 3]);
V_cubic(:,:,1) = flipud(fliplr(V8(:,:,2)));
V_cubic(:,:,2) = flipud(fliplr(V8(:,:,1)));
V_cubic


[Exx_c, Eyy_c, Ezz_c] = gradient(V_cubic, H, H, H);
E_cubic = {Exx_c, Eyy_c, Ezz_c};


%{
% Li passo el potencial.
% Camp als nodes de la malla fina que rodejen al meu punt:
E_is3 = cell(1,3);
for k = 1:3  % Cada componente del campo la promedio
    E_is = zeros(2,2,2);
    for pos = 1:8  % Itero matriu 2.2.2
        E_is(pos) = Efield{k}(ix(pos), iy(pos), iz(pos));
    end
    E_is3{k} = E_is;
end
%}
%end