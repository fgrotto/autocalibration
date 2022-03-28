function [best_H, best_focal_length, best_cost] = practical_autocal(keyframe1, keyframe2, P)
    % Extract the projective matrixes needed
    P1 = P(:,:,keyframe1);
    P2 = P(:,:,keyframe2);
    
    % Define the number of focal to try and the space
    num_focal_to_try = 100;
    focal_length_to_try = [];
    
    % Try focal length from 0.3 to 3.0 as suggested in the paper
    log_start = log(0.3);
    log_end = log(3.0);
    log_step = (log_end - log_start) / num_focal_to_try;
    for i = 1:num_focal_to_try
        focal_length_to_try(i) = exp(log_start+i*log_step);
    end
    
    % Iterate over each focal length to compute the related coste
    for i = size(focal_length_to_try)
        f = focal_length_to_try(i);
        
        % Update the homography from cameras given the focal length
        % Assuming that
        %             P1 = [Q | q]
        % then the transformation H, shown below,
        %             H = [  Q^-1  | -Q^-1 * q ] ,
        %                 [  0 0 0 |         1 ]
        % 
        % brings the camera P1 to the identity since
        %             [Q | q] * [  Q^-1  | -Q^-1 * q ] = [I | 0]
        %                       [  0 0 0 |         1 ]
  
        Q_inv = inv(P1(1:3,1:3));
        q = P1(:,4);
        
        % Transform both cameras such that P1 = [I|0].
        H_normalizer = [Q_inv, q; zeros(1,3), 1];
        P2H = P2*H_normalizer;
        
        % Compute plane at infinity
        % In the paper's equation 3, they use P2' = [Q2|q2], but instead take
        %
        %   P2' = [Q|q]
        %
        % for simplicity.
        Q = P2H(1:3,1:3);
        q = P2H(:,4);
        
        % the cheriality of the reconstruction is already accounted
        lambda = 1;
        
        % This assumes that K is the same for both camera, and centering
        % is already corrected
        K = [f, 0, 0; 0, f, 0; 0, 0, 1];
        K_inv = inv(K);
        
        % With K, q, and lambda, compute t2 for the second camera; note that the
        % transform by H made t1 = 0.
        t2 = K_inv * lambda * q;
        
        % First column is x; this will transform to (||x||, 0, 0).
        M(:,1) = t2;
        % Get a perpendicular vector of the previous one
        M(:,2) = [t2(3), t2(3), -t2(1) - t2(2)];
        % Get a perpendicular vector of the previous two
        M(:,3) = cross(t2, M(:,2));
        [Q_tmp,~] = qr(M);
        
        % This algorithm finds R (a rotation) such that
        % 
        % R * x == [ ||x|| 0 0 ].
        % 
        % If x is the vector, we can make a 3x3 matix [x * *], then take the QR
        % decomposition giving
        % M = [x * *] = QR
        % 
        % where, in this case, Q is a rotation and R is (confusingly) upper
        % triangular. Then, multiplying both sides with Q^-1 ( == Q^T)
        % 
        % Q^T M = R
        % 
        % Since the first column of R is of the form [y 0 0], and since rotations are
        % norm-preserving, y = ||x|| and therefore Q^T is the required matrix.
        R_star = Q_tmp';
        
        % Compute the W matrix based on equation (7), and extract w1..w3.
        W = R_star * K_inv * Q * K;
        w1 = W(1,:);
        w2 = W(2,:);
        w3 = W(3,:);
        
        % Compute v based on equation (9).
        v = (cross(w2, w3) / norm(w3) - w1) / norm(t2);
        
        % Build the updated homography
        H = [f, 0, 0, 0; 0, f, 0, 0; 0, 0, 1, 0; v, 1];
        H = H_normalizer * H;
        
        % Save H
        Hs{i} = H;
        
        % Compute the associated cost
        cost = 0;
        for j = 1:size(P,3)
           if j ~= keyframe1 && j ~= keyframe2 
               P_tmp = P(:,:,j)*H;
               cost = cost + cost_for_projective_matrix(f, P_tmp);
           end
        end
        
        % Save compute cost for the related H
        costs(i) = cost;
    end
    
    % Find the best cost, the best focal length and the best H
    best_cost = costs(1);
    for i = 1:size(focal_length_to_try)
        if (costs(i) <= best_cost)
          best_focal_length = focal_length_to_try(i);
          best_H = Hs{i};
          best_cost = costs(i);
        end
    end
end
