function [R,t] = relative_lin(q2, q1, K1, K2)
%RELATIVE_LIN Relative orientation

    q1 = K1\ensure_homogeneous(q1); 
    q2 = K2\ensure_homogeneous(q2); 
    % now q1 and q2 are homogeneous NIC
    
    E = eight_points(q2(1:2,:), q1(1:2,:)) ;
    
    [U,~,V] = svd(E);
    
    S1 = [0 1 0; -1 0 0; 0 0 0];
    R1 = [0 -1 0; 1 0 0; 0 0 1];
    
    for j = 1:4
        % 4 combinations
        S = (-1)^j * U*S1*U';
        
        if j==3   R1 = R1'; end
        
        t=[S(3,2) S(1,3) S(2,1)]';
        R = det(U*V')*U*R1*V';
        
        % 3D points must be in front of both cameras
        [z2, z1] =  icomb(q2(:,1), -R*q1(:,1),t );
        if z1>0 && z2>0  break,  end
        
    end
    
    
