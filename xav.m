%% FINAL S5

clear all;
close all;
clc;
% load('S5_pot1');
% Paràmetres físics: 
global epsilon q m kel Dt
epsilon = 8.85e-12;
q = 3*1.6*10^-19; 
m = 196.96657*(10^-3)/(6.022*10^23); % Au3+
kel = 1/(4*pi*epsilon);
Dt = 10^-4;

% Paràmetres principals del problema:
global Ls Ns V0 a_sp b_sp dim
Ls = 10; % Space is a cube of dimensions 2*Ls.
Ns = 15;  % each dimension discretized in Ns nodes
V0 = 3
a_sp = -Ls/2; b_sp=Ls/2;
dim = a_sp + (0:Ns)*(b_sp-a_sp)/Ns;
metode_aprox = 0; % 0: interp3, 1: ponderacio

% Number of time instants in wich the trajectory os calculated.
Tkin = 25000;
tfin = Dt*Tkin;
t = 0:Dt:tfin;

% Paràmetres que hem de modificar:
Nperiod = 45; % Amb aixo modifiquem la frequencia. més N és menys w
th = linspace(0,2*pi,2*Nperiod+1); Dth = th(2)-th(1);
th = th(1:end-1);
w = Dth/Dt
potencials = V0*sin(th);
%potencials = V0*sign(sin(th));
%potencials = -V0*sawtooth(th)
mod = repmat(potencials,[1,ceil(Tkin/(2*Nperiod))]);

% IF NECESSARY, RUN THE FUNCTION:
S5_pot


%% TRAJECTORIES

%{
Cell de matrius de condicions inicials:
M = [x0, y0, z0; x0', y0', z0'; x0'', y0'', z0''];
CI = {M1, M2... M_Npart};
%}
Npart = 10;
CI = cell(Npart,1);
CI{1} = [1e-3*rand 1e-3*rand 1e-14*rand; -5*rand -5*rand -1e-14*rand; 0 0 0];
CI{2} = [-1e-3*rand -1e-3*rand -1e-14*rand; 5*rand 5*rand 1e-14*rand; 0 0 0];
CI{3} = [-1e-3*rand 1e-3*rand -1e-14*rand; 5*rand -5*rand -1e-14*rand; 0 0 0];
CI{4} = [-1e-3*rand -1e-3*rand -1e-14*rand; -5*rand 5*rand -1e-14*rand; 0 0 0];
CI{5} = [-1e-3*rand -1e-3*rand 1e-14*rand; -5*rand -5*rand 1e-14*rand; 0 0 0];
CI{6} = [1e-3*rand 1e-3*rand 1e-14*rand; -5*rand -5*rand -1e-14*rand; 0 0 0];
CI{7} = [-1e-3*rand -1e-3*rand -1e-14*rand; 5*rand 5*rand 1e-14*rand; 0 0 0];
CI{8} = [-1e-3*rand 1e-3*rand -1e-14*rand; 5*rand -5*rand -1e-14*rand; 0 0 0];
CI{9} = [-1e-3*rand -1e-3*rand -1e-14*rand; -5*rand 5*rand -1e-14*rand; 0 0 0];
CI{10} = [-1e-3*rand -1e-3*rand 1e-14*rand; -5*rand -5*rand 1e-14*rand; 0 0 0];




% Omplo les cells perque tinguin ions
POS = cell(Npart, 1);
X = zeros(Tkin-1,3);
for n = 1:Npart
    POS{n} = X;
end
VEL = POS; ACC = POS; F = POS;

for n = 1:Npart
    POS{n}(1,:) = CI{n}(1,:);
    VEL{n}(1,:) = CI{n}(2,:);
    ACC{n}(1,:) = CI{n}(3,:);
end

Vt_rep = ones(1,Tkin-1); %REP
for i = 1:Tkin-1
    indraw = zeros(1,3);
    for k = 1:3 % get index for Efield
        indraw(k) = (Ns)*(POS{n}(i,k)-a_sp)/(b_sp-a_sp)+1;
        if indraw(k) > length(dim) - 2 || indraw(k) < 2
            disp('!!')
            disp('LA PARTICULA MARXA DE LA REGIO')
            disp('!!')
        end
    end
    for n = 1:Npart
        [E, V_per_representar] = linear_E(POS{n}(i,:), dim, indraw, mod(i), V, metode_aprox);
        Vt_rep(i) = V_per_representar; %REP
        
        % Repulsion between particles:
        rn = POS{n}(i,:);
        for p = [1:n-1, n+1:Npart] % evito la propia
            rp = POS{p}(i,:); R = rp-rn;
            F_pn = kel*(q^2).*R./(norm(R)^3);
            F{n}(i,:) = F{n}(i,:) + F_pn;
        end
        
        for k = 1:3
            ACC{n}(i,k) = (q/m)*E(k) + F{n}(i,k);
            VEL{n}(i+1,k) = VEL{n}(i,k) + ACC{n}(i,k)*Dt;
            POS{n}(i+1,k) = POS{n}(i,k) + VEL{n}(i,k)*Dt + 0.5*ACC{n}(i,k)*(Dt^2);
        end
    end
end


[Mvertex, Mtopol, Nt_lat] = object_ion_trap(Ntop, Nlat, ztop, Ro_lat, h, Ls, true); hold on;
X = POS{1};
plot3(X(1:i,1),X(1:i,2),X(1:i,3)); hold on;
X = POS{2};
plot3(X(1:i,1),X(1:i,2),X(1:i,3));
X = POS{3};
plot3(X(:,1),X(:,2),X(:,3));
hold off;

figure;
X = POS{1};
plot3(X(:,1),X(:,2),X(:,3)); hold on;
title('Trajectory')
X = POS{2};
plot3(X(:,1),X(:,2),X(:,3));
X = POS{3};
plot3(X(:,1),X(:,2),X(:,3));
xlabel('x'); ylabel('y'); zlabel('z');
hold off;

figure;
sgtitle('Els 4 imprescindibles');
subplot(2,2,1);
plot(t(1:end-2), Vt_rep); hold on;
title('Potencial al centre');
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