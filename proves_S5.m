%% WORKS S3.5

clear all;
close all;
clc;

epsilon = 8.85e-12; % Si la canvio aqui la he de canviar tmb a potential_in_space
Ls = 10;  % space is a cube of side 2*Ls
Ns = 15;  % each dimension discretized in Ns nodes

Tkin = 5000;    % Number of time instants in wich the trajectory os calculated.

% NOVES VARIABLES:
%V0=0.1;
V0 = 10;
a_sp = -Ls; b_sp=Ls; j = 0:Ns;

Npart = 2;

% Per les dues maneres de plantejar el problema, tenim:
aprox = true; % aproximo el camp electric per diferents metodes
%--------------------------------------------------------------------------
% METODO DE LOS 1, -1 FUNCIONA
Dt = 10^-4;
q = 1.6*10^-19; 
m = 196.96657*(10^-3)/(6.022*10^23); % Au3+
kel = 1/(4*pi*epsilon);

Nperiod = 25; 
th = linspace(0,2*pi,2*Nperiod+1); Dth = th(2)-th(1);
th = th(1:end-1);
w = Dth/Dt
potencials = V0*sin(th);
mod = repmat(potencials,[1,ceil(Tkin/(2*Nperiod))]);

%potencials = V0*sign(sin(th));
%potencials = -V0*sawtooth(th)

tfin = Dt*Tkin;
t = 0:Dt:tfin;


% DISCRETIZATION OF ALL SPACE:
% Malla gruixuda
dim = a_sp + j*(b_sp-a_sp)/Ns; H = dim(2)-dim(1); % espaciado entre puntos
[xx,yy,zz] = meshgrid(dim, dim, dim); msh = {xx, yy, zz};

% HERE WE DISCRETIZE DESIRED GEOMETRY OF THE PLATES:
Ntop = 7; Nlat = 4;
ztop = 3;
Ro_lat = 7;
h = 1;
[Mvertex, Mtopol, Nt_tops] = object_ion_trap(Ntop, Nlat, ztop, Ro_lat, h, Ls, false);

%--------------------------------------------------------------------------
% AQUI ES ON CANVIA EL METODE: NOMES CALCULO UN POTENCIAL.
% POTENTIAL AT THE PLATES:
Vr = potential_in_space(1, Mvertex, Mtopol, Nt_tops, msh);
V = Vr{1};


%% TRAJECTORIES

% Omplo les cells perque tinguin ions
POS = cell(Npart, 1);
X = zeros(Tkin-1,3);
for n = 1:Npart
    POS{n} = X;
end
VEL = POS; ACC = POS; F = POS;

% Initial conditions:
% POS{1}(1,1) = 0.3; POS{2}(1,1) = 0.1;
% POS{1}(1,1) = (rand+0.5); POS{1}(1,2) = (4*rand-2); POS{1}(1,3) = (rand+0.5);
% VEL{1}(1,1) = (-6); VEL{1}(1,1) = (-5); VEL{1}(1,1) = (1);
% POS{2}(1,1) = (2*rand-3); POS{2}(1,2) = (-3*rand+1); POS{2}(1,3) = -rand+0.05;
% VEL{2}(1,1) = (5); VEL{2}(1,2) = (10); VEL{2}(1,3) = (2);

POS{1}(1,1) = (rand); POS{1}(1,2) = (-rand); POS{1}(1,3) = (10^-12*rand);
VEL{1}(1,1) = (-100*rand); VEL{1}(1,2) = (90*rand); VEL{1}(1,3) = (-50*rand);
POS{2}(1,1) = (-0.3); POS{2}(1,2) = (-0.3); POS{2}(1,3) = (10^-12*-rand);
VEL{2}(1,1) = (25*rand); VEL{2}(1,2) = (15*-rand*10-4); VEL{2}(1,3) = (rand*6-4);


Vt_rep = ones(1,Tkin-1); %REP
E_rep = zeros(Tkin-1,3); %REP
Eap = E_rep; % REP
for i = 1:Tkin-1
    indraw = zeros(1,3);
    for k = 1:3 % get index for Efield
        indraw(k) = (Ns)*(POS{n}(i,k)-a_sp)/(b_sp-a_sp)+1;
        if indraw(k) > length(dim) || indraw(k) < 1
            disp('!!')
            disp('LA PARTICULA MARXA DE LA REGIO')
            disp('!!')
        end
    end
    for n = 1:Npart
        [Eap(i,:), V_per_representar] = linear_E(POS{n}(i,:), dim, indraw, mod(i), V);
        Vt_rep(i) = V_per_representar; %REP
        
        % Repulsion between particles:
        rn = POS{n}(i,:);
        for p = [1:n-1, n+1:Npart] % evito la propia
            rp = POS{p}(i,:); R = rp-rn;
            F_pn = kel*(q^2).*R./(norm(R)^3);
            F{n}(i,:) = F{n}(i,:) + F_pn;
        end
        
        for k = 1:3
            ACC{n}(i,k) = (q/m)*Eap(i,k) + F{n}(i,k);
            VEL{n}(i+1,k) = VEL{n}(i,k) + ACC{n}(i,k)*Dt;
            POS{n}(i+1,k) = POS{n}(i,k) + VEL{n}(i,k)*Dt + 0.5*ACC{n}(i,k)*(Dt^2);
        end
    end
end


[Mvertex, Mtopol, Nt_lat] = object_ion_trap(Ntop, Nlat, ztop, Ro_lat, h, Ls, true); hold on;
X = POS{1};
plot3(X(1:i,1),X(1:i,2),X(1:i,3)); hold on;
% X = POS{2};
% plot3(X(:,1),X(:,2),X(:,3));
hold off;

figure;
X = POS{1};
plot3(X(:,1),X(:,2),X(:,3)); hold on;
title('Trajectory')
% X = POS{2};
% plot3(X(:,1),X(:,2),X(:,3));
xlabel('x'); ylabel('y'); zlabel('z');
hold off;

figure;
sgtitle('Trajectory');
subplot(2,2,1);
plot(t(1:end-2), Vt_rep); hold on;
title('Potential at the center');
subplot(2,2,2);
plot(t(1:end-1), VEL{1}(:,1)); hold on;
title('vx(t)');
subplot(2,2,3);
plot(t(1:end-2), ACC{1}(:,1)); hold on;
title('ax(t)');
subplot(2,2,4);
plot(t(1:end-1), POS{1}(:,1)); hold on;
title('x(t)');
hold off;

figure;
sgtitle('Els altres imprescindibles');
subplot(3,1,1);
plot(t(1:end-2), Vt_rep); hold on;
title('Potencial al centre');
subplot(3,1,2);
plot(t(1:end-2), Eap(:,1)); hold on;
title('Exap(t)');
% subplot(3,1,3);
% plot(t(1:end-2), E_orig(:,1)); hold on;
% title('Ex(t) original');
hold off;

figure;
for n = 1:Npart
    X = POS{n};
    subplot(3,1,1);
    plot(t(1:end-1), X(:,1));
    title('x'); hold on;
    subplot(3,1,2);
    plot(t(1:end-1), X(:,2));
    title('y'); hold on;
    subplot(3,1,3);
    plot(t(1:end-1), X(:,3));
    title('z'); hold on;
end
hold off;

% figure;
% col = {'r', 'k'};
% for i = 1:Tkin-1
%     for n = 1:Npart
%         X = POS{n};
%         plot3(X(i,1),X(i,2),X(i,3), '-o','Color', col{n}); hold on;
%         axis([-1 1.5 -2 2 -0.7 0.8]);
%     end
%     pause(10^-16);
% end
% hold off;

V_ini = ipermute(V, [3 1 2]); % Faig un tall del potencial (el vull a x=cnt)
V_ini = V_ini(:,:,ceil(Ns/2)); % tall a x = cnt
[xx, yy] = meshgrid(dim);
figure;
surf(xx,yy,V_ini);
colormap jet; shading interp; colorbar;
axis([-Ls Ls -Ls Ls -V0/2 V0/2]);
xlabel('x'); ylabel('y'); zlabel('V(t=0, y=0)');
title('V(t=0, y=0) surface');
hold on;
plot3(X(:,1), X(:,2), X(:,3)+V_ini(3,3), 'k');
hold off;