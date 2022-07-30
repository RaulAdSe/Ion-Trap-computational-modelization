function [E_lin, V_per_representar] = linear_E(X, dim, indraw, mod, Vi, metode)
               
mi = zeros(3,2);
for k = 1:3
    if floor(indraw(k)) == ceil(indraw(k))
        indraw(k) = indraw(k) + 0.1; % No passa mai pero per si de cas
    end
    mi(k, 1) = floor(indraw(k))-1;
    mi(k, 2) = ceil(indraw(k))+1;
end
[ix, iy, iz] = meshgrid([mi(1,1):mi(1,2)], [mi(2,1):mi(2,2)], [mi(3,1):mi(3,2)]);
[mx, my, mz] = meshgrid([dim(mi(1,1):mi(1,2))], [dim(mi(2,1):mi(2,2))], [dim(mi(3,1):mi(3,2))]);
mmsh = {mx, my, mz};
H = dim(2) - dim(1);
V_cubic = zeros(4,4,4);
for k = 1:64
    V_cubic(k) = Vi(ix(k), iy(k), iz(k));
end
[Exx_c, Eyy_c, Ezz_c] = gradient(-mod*V_cubic, H); E_c = {Exx_c, Eyy_c, Ezz_c};

if metode == 0
    E_lin = zeros(1,3);
    E_lin(1) = interp3(mx,my,mz,Exx_c,X(1),X(2),X(3));
    E_lin(2) = interp3(mx,my,mz,Eyy_c,X(1),X(2),X(3));
    E_lin(3) = interp3(mx,my,mz,Ezz_c,X(1),X(2),X(3));
    
elseif metode == 1     
    sums = zeros(4,4,4);
    for r = 1:3
        mr = mmsh{r} - X(r);
        sums = sums + mr.^2;
    end
    sums = fliplr(sqrt(sums));

    punt_coincideix = false;
    for k = 1:64  % Si la particula esta just a un node, no cal ponderar
        if sums(k) < 10^-14 % aixo per mi ja es suficientment aprop.
            punt_coincideix = true;
            value_coinc = k;
        end
    end

    if punt_coincideix == true
        E_pond = [E_c{1}(value_coinc), E_c{2}(value_coinc), E_c{3}(value_coinc)];

    else % Si no, he de fer la mitjana ponderada.
        pesos = (1./sums).^2; % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        E_pond = zeros(1,3);
        for k = 1:3
            sum_pesos = sum(pesos,'all');
            E_pond(k) = (1/sum_pesos)*sum(E_c{k}.*pesos,'all'); % Mateixa sumd per tots.
            if sum_pesos < 10^(-15)
                E_pond(k) = 0;
            end
        end
    end
    E_lin = E_pond;
end

V_per_representar = mod*V_cubic(1,1,1);
end
