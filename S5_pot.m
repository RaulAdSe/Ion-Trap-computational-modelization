%% CALCUL DEL POTENCIAL A TOT L'ESPAI

close all;
global Ls Ns a_sp b_sp 

% DISCRETIZATION OF ALL SPACE:
% Malla gruixuda única
dim = a_sp + (0:Ns)*(b_sp-a_sp)/Ns;
[xx,yy,zz] = meshgrid(dim, dim, dim); msh = {xx, yy, zz};

% HERE WE DISCRETIZE DESIRED GEOMETRY OF THE PLATES:
Ntop = 7; Nlat = 4;
ztop = 3;
Ro_lat = 7;
h = 1;
[Mvertex, Mtopol, Nt_tops] = object_ion_trap(Ntop, Nlat, ztop, Ro_lat, h, Ls, false);

Vr = potential_in_space(1, Mvertex, Mtopol, Nt_tops, msh);
V = Vr{1};

save('S5_pot');