function [Mvertex, Mtopol, Nt_tops] = object_ion_trap(Ntop, Nlat, ztop, Ro_lat, h, Ls, rep)
% OBJECT_ION_TRAP Generates the object of the geometry of Paul ion trap.
% Inputs:
    % Ntop: number of concentric circles for the top geometry
    % Nlat: number of concentric circles for the lateral geometry
    % ztop: initial height for top geometry
    % Ro_lat: initial radius for lateral geometry
    % h: step
    % Ls: space region
    % rep: boolean
% Outputs:
    % Mvertex: Matrix containing the vertex of the triangularization of the geometry.
    % Mtopol: Matrix assigning the vertex to a triangle.
    % Ntpos: Number of triangles of one the "parts". It will be subject to
    % the potential Viter/2. Ntneg = Nt-Ntpos will be the other part,
    % subject to -Viter/2.

%% Hiperboloide 2 fulles

N = Ntop; % N - 1 concentric circles (first is point r = 0)
rconc = zeros(N,1);
for i = 2:N
    rconc(i) = rconc(i-1) + h;
end
% arch of circumference = h, equispaced points for each circle: phi = h/r
steps = round(2*pi*rconc/h); 
htrue = (2*pi*rconc)./steps; htrue(1) = 0;

% Geometry points
xx = 0; yy = 0; zz = ztop;
for i = 2:N
    n = (1: 1: steps(i));
    xx = [xx rconc(i)*cos(n*htrue(i)/rconc(i))];
    yy = [yy rconc(i)*sin(n*htrue(i)/rconc(i))];
    z = sqrt(0.5)*sqrt(2*ztop^2 + rconc(i)^2);
    zr = repmat(z, 1, steps(i));
    zz = [zz zr];
end

% Triangulation
obj.vertex = [xx(:), yy(:), zz(:)]'; % top
obj.vertex2 = [xx(:), yy(:), -zz(:)]'; % bottom
Mvertex = [obj.vertex obj.vertex2];

obj.topol = delaunay(xx(:), yy(:))';
Mtopol = [obj.topol obj.topol + max(obj.topol, [], 'all')];
Nt_tops = length(Mtopol); % Triangles built in the tops.

% Representation
if rep == true 
    vx=zeros(size(Mtopol,2),3); vy = vx; vz = vx;
    figure;
    for i=1:size(Mtopol,2)
        vx(i,:)=Mvertex(1,Mtopol(:,i));
        vy(i,:)=Mvertex(2,Mtopol(:,i));
        vz(i,:)=Mvertex(3,Mtopol(:,i));
        fill3(vx(i,:), vy(i,:), vz(i,:), 'g');
        hold on
    end
end

%% Hiperboloide 1 fulla

ss = round(2*pi*Ro_lat/h); % each circle has the same amount of samples
hrr = 2*pi*Ro_lat/ss; % each circle has its step
N = Nlat; 

rconc = zeros(N,1);
% First circle
rconc(1) = Ro_lat;
n = (1: 1: ss);
xx = rconc(1)*cos(n*hrr/rconc(1));
yy = rconc(1)*sin(n*hrr/rconc(1));
z = sqrt(0.5)*sqrt(rconc(1)^2 - Ro_lat^2);
zr = repmat(z, 1, ss);
zz = zr;

% function fun(R) to find new proper radii
tol = 1e-4; itmax = 15; a = 0.001;
for i = 2:N
    fun = @(R)(R.^4 + 4/3*rconc(i - 1)*R.^3 + 4/9*(3*rconc(i - 1)^2 - 3*hrr^2 - 2*Ro_lat^2)*R.^2 - 8/9*rconc(i - 1)*hrr^2*R + 4/9*hrr^2*(1 - 2*rconc(i - 1)^2 + 2*Ro_lat^2));
    [xk,~] = newton(a, tol, itmax, fun);
    a = a^a;
    rconc(i) = xk(end) + rconc(i - 1);
    hrr = 2*pi*rconc(i)/ss;
    
    n = (1: 1: ss);
    xx = [xx rconc(i)*cos(n*hrr/rconc(i))];
    yy = [yy rconc(i)*sin(n*hrr/rconc(i))];
    z = sqrt(0.5)*sqrt(rconc(i)^2 - Ro_lat^2);
    zr = repmat(z, 1, ss);
    zz = [zz zr];
end

% Triangulation
obj.vertex1 = [xx(:), yy(:), zz(:)]';
obj.vertexn = [xx(:), yy(:), -zz(:)]';
Mv_lat = [obj.vertex1 obj.vertexn]; % new vertex matrix

obj.topol_raw = delaunay(xx(:), yy(:))'; % 1:ss have z = 0
obj.topol = []; % we get rid of triangles build upon z = 0 plane
for i = 1:size(obj.topol_raw, 2)
    if obj.topol_raw(1, i) > ss || obj.topol_raw(2, i) > ss || obj.topol_raw(3, i) > ss
        obj.topol = [obj.topol obj.topol_raw(:,i)];
    end
end
obj.topol2 = obj.topol +  max(obj.topol, [], 'all'); % conditioned topology indexes

Mt_lat = [obj.topol obj.topol2]; % new topology matrix

% Representation
if rep == true 
    vx=zeros(size(Mt_lat,2),3); vy = vx; vz = vx; 
    for i=1:size(Mt_lat,2)
        vx(i,:) = Mv_lat(1,Mt_lat(:,i));
        vy(i,:) = Mv_lat(2,Mt_lat(:,i));
        vz(i,:) = Mv_lat(3,Mt_lat(:,i));
        fill3(vx(i,:), vy(i,:), vz(i,:), 'g');
        hold on
   end
    axis([-Ls Ls -Ls Ls]);
    xlabel('x'); ylabel('y'); zlabel('z');
    title('Geometry Plot');
end
%hold off;
hold on;

% To merge two different topology matrixes, one's indexes have to be conditioned
Mt_lat = Mt_lat + max(Mtopol, [], 'all'); 
Mvertex = [Mvertex Mv_lat];
Mtopol = [Mtopol Mt_lat];
