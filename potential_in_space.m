function Vr = potential_in_space(Viter, Mvertex, Mtopol, Ntpos, msh)
%POTENTIAL_IN_SPACE It calculates the potential in desired points of space for a given discretization of a metallic geometry subject to a certain potential.
%{
Regardless the geometry, it must consist of two "parts" or components, one
of them subject to a given potential (Viter/2), and the other to the opposite
potential (-Viter/2), in order to generate a saddle point that can confine ions in
its interior.
This function assumes both "parts" are included in the geometry
matrices (topology and vertex), the first consisting of Ntpos triangles, and 
the second of Ntneg triangles.
%}
% Inputs:
    % Viter: time-varying potential in the plates' geometry.
    % Mvertex: Matrix containing the vertex of the triangularization of the geometry.
    % Mtopol: Matrix assigning the vertex to a triangle.
    % Ntpos: Number of triangles of one the "parts". It will be subject to
    % the potential Viter/2. Ntneg = Nt-Ntpos will be the other part,
    % subject to -Viter/2.
    % msh: Cell {xx,yy,zz} containing the 3D discretization of all space.

% Outputs:
    % Vr: a cell of size 1xT with the potential at the mesh (msh).

epsilon = 8.85e-12;
    
% Coses que m'estalvio passar com a variable:
T = length(Viter); % how many potentials I must calculate. CC increases.
Nt = size(Mtopol, 2); % total number of triangles
Ntneg = Nt-Ntpos;
Ns = size(msh{1},1);

% Position of all vertexes i.Nmf
v1 = Mvertex(:,Mtopol(1,:));
v2 = Mvertex(:,Mtopol(2,:));
v3 = Mvertex(:,Mtopol(3,:));
obj.cent =(v1+v2+v3)/3; % Center of the triangle

% Unitary surface vector for all triangles:
r12 = v2-v1;
r13 = v3-v1;
c = cross(r12,r13); % Already by columns
obj.norm_S = sqrt(sum(c.^2));   % Makes the normal o
normalmat = repmat(obj.norm_S,3,1);
obj.unit_S = c./normalmat;

% Finally we can construct matrix Z:
Z = zeros(Nt); % Int me dona la contribució de cada triangle de càrrega (fun base) a un trosset d'espai determinat Rf. Me dona files de Z
for n=1:Nt
    Rf = obj.cent(:,n); %Single point to compute integral vs centers (all)
    Z(n, :) = int_S_1divR(Rf, v1, v2, v3, obj.unit_S, obj.cent)/(4*pi*epsilon*obj.norm_S(n)/2);
end

Vr = cell(T,1);
for i=1:T % At every instant of time
    Vt = Viter(i);
    Vo = [ones(Ntpos,1)*(Vt/2); ones(Ntneg,1)*(-Vt/2)]; % top geometry and lateral geometry potentials
    Q = Z\Vo;
    
    %Potential at every point of space
    V = zeros(Ns,Ns,Ns);
    for j = 1:Ns^3
        Rf=[msh{1}(j),msh{2}(j),msh{3}(j)]';
        V(j)=int_S_1divR(Rf, v1, v2, v3, obj.unit_S, obj.cent)*Q*2/(4*pi*epsilon*obj.norm_S(n)); %int_S_1divR*Q returns a single value
    end
    Vr{i} = V;
end





