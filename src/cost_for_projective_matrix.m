function [cost] = cost_for_projective_matrix(focal_length, P)
    K = P(1:3,1:3);
    
    % Ensure that the K diagonal is positive
    if (K(3,3) < 0)
        K = -K;
    end
    if (K(2,2) < 0) 
        S = [1,  0,  0;
             0, -1,  0;
             0,  0,  1];
        K = K * S;
    end
    if (K(1,1) < 0)
        S = [-1, 0, 0;
              0, 1, 0;
              0, 0, 1];
        K = K * S;
    end
    
    % Scale to have K(3,3) = 1
    K = K/K(3,3);
    
    % Extract parameters from intrinsic camera parameter matrix
    fx = abs(K(1,1));
    fy = abs(K(2,2));
    avg_focal_length = 0.5 * (fx + fy);
    skew = abs(K(1,2));
    cx = abs(K(1,3));
    cy = abs(K(2,3));
    
    % Penalize deviation from the guessed focal length.
    if max(focal_length, avg_focal_length) == 0
        focal_length_error_cost = 0;
    else
        M = max(focal_length, avg_focal_length);
        m = min(focal_length, avg_focal_length);
        focal_length_error_cost = 1 - (M-m) / M;
    end
    
    % Penalize deviation from equal aspect ratio.
    if max(fx, fy) == 0
        fraction_off_equal_aspect_ratio_cost = 0;
    else
        M = max(fx, fy);
        m = min(fx, fy);
        fraction_off_equal_aspect_ratio_cost = 1 - (M-m) / M;
    end
    
    % Penalize non-zero center of projection.
    principal_point_cost = cx+cy;
    
    cost = 1. / 0.05   * focal_length_error_cost + 1. / 0.001  * fraction_off_equal_aspect_ratio_cost+ 1. / 0.0005 * skew + 1. / 0.01   * principal_point_cost;
    
    % Rescale actual cost
    cost = cost / 100; 
    cost = cost * cost; 
end