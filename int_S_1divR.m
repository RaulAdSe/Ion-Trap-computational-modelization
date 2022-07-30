% Computation of the 1/R singular Kernel contribution of a set of triangles (r_1,r_2,r_3) in rf: 
%
%  Output: tmp_1divR; intS(1/R)
%
%  Input: rf: field point [3x1]
%         r_1 : first vertex of triangles [3xNt]
%         r_2 : second vertex of triangles [3xNt]
%         r_3 : third vertex of triangles [3xNt]
%         un_S: normal vectors associated to the set of triangles <r1-r2-r3> [3xNt]
%         < IMPORTANT: r1-r2-r3 must turn in the same sense as un_S  >
%         cent_S: set of centroids of the triangles <r1-r2-r3> [3xNt]
%
% by Ed. Ubeda, March 2008

function  tmp_1divR = int_S_1divR( rf , r_1 , r_2 , r_3 , un_S , cent_S);

N = size(r_1,2);

tmp_1divR = zeros(1,N);

% Tangential vectors
l_vec_12 = unitary( r_2 - r_1 );
l_vec_23 = unitary( r_3 - r_2 );
l_vec_31 = unitary( r_1 - r_3 );

% Normal-contour vectors
un_c_vec_12 = cross( l_vec_12 , un_S );
un_c_vec_23 = cross( l_vec_23 , un_S );
un_c_vec_31 = cross( l_vec_31 , un_S );

% Normal signed-distance from rf to the surface of the triangles
D_S = sum( ( rf*ones(1,N)  - cent_S ) .* un_S );

for s=1:3,
    
    if (s==1),
        %% 12
        v_plus = r_2;
        v_min = r_1;
        un_c = un_c_vec_12;
        l_vec = l_vec_12;
    elseif (s==2),
        %% 23
        v_plus = r_3;
        v_min = r_2;        
        un_c = un_c_vec_23;
        l_vec = l_vec_23;
    elseif (s==3),
        %% 31
        v_plus = r_1;
        v_min = r_3;                
        un_c = un_c_vec_31;
        l_vec = l_vec_31;
    end;  %%% if (s==1),
    

    D_vec_plus = ( ones(3,1) * sum( ( v_plus - rf*ones(1,N) ) .* un_S ) ) .* un_S;
    rho_vec_plus = ( v_plus - rf*ones(1,N) ) - D_vec_plus;
    
    D_vec_min = ( ones(3,1) * sum( ( v_min - rf*ones(1,N) ) .* un_S ) ) .* un_S;    
    rho_vec_min  = ( v_min  - rf*ones(1,N) ) - D_vec_min;     


    l_plus = sum( rho_vec_plus .* l_vec );
    l_min = sum( rho_vec_min .* l_vec );        
    
    Po = abs( sum( rho_vec_plus .* un_c ) );
    
    Po_other = abs( sum( rho_vec_min .* un_c ) );
    
    i_Po_ok = find( Po > eps );
    Po_vec = zeros(3,N);
    Po_vec(:,i_Po_ok) = ( rho_vec_plus(:,i_Po_ok) - ( ones(3,1)*l_plus(i_Po_ok) ) .* l_vec(:,i_Po_ok) ) ./ ( ones(3,1) * Po(i_Po_ok) );
    
    P_plus = sqrt( Po.^2 + l_plus.^2 );
    P_min = sqrt( Po.^2 + l_min.^2 );
    
    Ro = sqrt( Po.^2 + D_S.^2 );
    R_plus = sqrt( P_plus.^2 + D_S.^2);
    R_min = sqrt( P_min.^2 + D_S.^2);    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Computation of the scalar term: 1/R %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% term_1 : Po*log() term _plus, _min    
    term_1 = zeros(1,N);
    
    term_1_plus_log = R_plus + l_plus;    
    i_ok_plus_log = find( term_1_plus_log > eps );
    term_1(i_ok_plus_log) = + Po(i_ok_plus_log) .* log( term_1_plus_log( i_ok_plus_log ) );

    term_1_min_log = R_min + l_min;        
    i_ok_min_log = find( term_1_min_log > eps );    
    term_1(i_ok_min_log) = term_1(i_ok_min_log) - Po(i_ok_min_log) .* log( term_1_min_log( i_ok_min_log ) );        
    
    %% term_2 : atan() term _plus, _min
    term_2 = zeros(1,N);
    
    term_2_plus_a = Po .* l_plus;
    term_2_plus_b = Ro.^2 + abs(D_S).*R_plus;
    i_ok_plus_ab = find( term_2_plus_b > eps );
    term_2(i_ok_plus_ab) = - abs( D_S(i_ok_plus_ab) ) .* atan( term_2_plus_a(i_ok_plus_ab)./term_2_plus_b(i_ok_plus_ab) );
    
    term_2_min_a = Po .* l_min;
    term_2_min_b = Ro.^2 + abs(D_S).*R_min;    
    i_ok_min_ab = find( term_2_min_b > eps );    
    term_2(i_ok_min_ab) = term_2(i_ok_min_ab) + abs( D_S(i_ok_min_ab) ) .* atan( term_2_min_a(i_ok_min_ab)./term_2_min_b(i_ok_min_ab) );    
    
    tmp_1divR = tmp_1divR + sum( Po_vec .* un_c ) .* ( term_1 + term_2 );
    
end; %%% for s=1:3,