function [K_cipolla] = compute_mc(Fs, K0)
    % Calculate intrinsic camera parameters
    K_vector = lsqnonlin(@(X) cost_medonca_cipolla(Fs, X),[K0(1,:) K0(2,2:3)],[],[], optimoptions('lsqnonlin','Display','off','Algorithm','levenberg-marquardt','TolX', 1e-12));
    
    % intrinsic camera parameter in matrix form
    K_cipolla = [K_vector(1) K_vector(2) K_vector(3); 
                0 K_vector(4) K_vector(5); 
                0 0 1];
end