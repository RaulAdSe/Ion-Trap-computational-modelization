%% FINAL S5

clear all;
close all;
clc;

% Paràmetres físics: 
global epsilon q m kel Dt
epsilon = 8.85e-12;
q = 1.6*10^-19; 
m = 196.96657*(10^-3)/(6.022*10^23); % Au3+
kel = 1/(4*pi*epsilon);
Dt = 10^-4;

% Paràmetres principals del problema:
global Ls Ns V0 a_sp b_sp dim
Ls = 10; % Space is a cube of dimensions 2*Ls.
Ns = 15;  % each dimension discretized in Ns nodes
V0 = 10;
a_sp = -Ls/2; b_sp=Ls/2;
dim = a_sp + (0:Ns)*(b_sp-a_sp)/Ns;
metode_aprox = 0; % 0: interp3, 1: ponderacio

% Number of time instants in wich the trajectory os calculated.
Tkin = 10000;  
tfin = Dt*Tkin;
t = 0:Dt:tfin;

% Paràmetres que hem de modificar:
Nperiod = 8; % Amb aixo modifiquem la frequencia
th = linspace(0,2*pi,2*Nperiod+1); Dth = th(2)-th(1);
th = th(1:end-1);
w = Dth/Dt
potencials = V0*sin(th);
%potencials = V0*sign(sin(th));
%potencials = -V0*sawtooth(th)
mod = repmat(potencials,[1,ceil(Tkin/(2*Nperiod))]);

% IF NECESSARY, RUN THE FUNCTION:
S5_pot
load('S5_pot');

%% TRAJECTORIES

%{
Cell de matrius de condicions inicials:
M = [x0, y0, z0; x0', y0', z0'; x0'', y0'', z0''];
CI = {M1, M2... M_Npart};
%}
Npart = 10;
CI = cell(Npart,1);
% CI{1} = [rand -rand rand; -170 30 5e-12 ; 0 0 0];
% CI{1} = [rand rand 0.1*rand; 200 -60 10e-11*rand; 0 0 0];
% CI{2} = [-0.3*rand -0.1*rand -0.1*rand; -5*rand 5*rand -20e-14*rand; 0 0 0];
% CI{3} = [0.1*rand -0.2*rand 0.05*rand; -4*rand 7*rand -10e-13*rand; 0 0 0];
% CI{4} = [0.2*rand 0.15*rand 0.1*rand; 8*rand -6*rand -10e-14*rand; 0 0 0];
% CI{5} = [-0.13*rand -0.21*rand -0.1*rand; 15*rand -8*rand -20e-14*rand; 0 0 0];
CI{1} = [0.0109697464523194	0.00635913709751057	0.0404579995857626; 1.34511873619949	1.46326470735269	7.63504640848813e-13; 0 0 0];
CI{2} = [-0.0627896379614169	-0.0771980385554245	-0.0932853570278820; -4.86370427001507	0.960141747138875	-2.77748405658311e-14; 0 0 0];
CI{3} = [0.0696266337082995	-0.00938200267748656	0.0262702201929668; -2.12137687357145	6.02797867975333	-4.84853333552102e-13; 0 0 0];
CI{4} = [0.0393456361215266	0.0671431139674026	0.0741257943454207; 4.16041973912309	-2.08627602766515	-1.49997253831683e-14; 0 0 0];
CI{5} = [-0.0586092067231462	-0.0262145317727807	-0.00444540922782385; 11.3239990084677	-1.94228286256770	-8.84804626003887e-14; 0 0 0];
CI{6} = [0.0687796085120107	-0.0359228210401861	0.0368170037150601; -1.57882990111505	10.2512380045197	-7.04047430334266e-13; 0 0 0];
CI{7} = [0.0442305413383371	0.00195776235533187	0.0330857880214071; 3.39447597466510	0.810811270296196	1.97053798095456e-14; 0 0 0];
CI{8} = [-0.0821721184961310	-0.0429921409383266	-0.0887770954256354; 5.86774493191745	-6.15291509910637	-7.93583034027233e-11; 0 0 0];
CI{9} = [0.0808514095887345	-0.0755077099007084	0.0188697772417551; -0.864075663845577	6.32325774373531	9.49303911849797e-14; 0 0 0];
CI{10} = [0.0327565434075205	0.0671264370451740	0.0438644982586956; -4.16750297794488	4.61312551457769	1.67253545494722e-10; 0 0 0];


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
        if indraw(k) > length(dim) || indraw(k) < 1
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
hold on
for i=1:Npart
X = POS{i};
plot3(X(:,1),X(:,2),X(:,3),'linewidth',1); hold on;
end
hold off;

figure;
for i=1:Npart
    X = POS{i};
plot3(X(:,1),X(:,2),X(:,3),'linewidth',1); hold on;
end
title('Trajectory')
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
    plot(t(1:end-1), X(:,1),'linewidth',1.5);
    title('x'); hold on;
    subplot(3,1,2);
    plot(t(1:end-1), X(:,2),'linewidth',1.5);
    title('y'); hold on;
    subplot(3,1,3);
    plot(t(1:end-1), X(:,3),'linewidth',1.5);
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