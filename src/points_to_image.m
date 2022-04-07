function [uv, vv, K, R, t] = points_to_image(filename, Points) 
    % Read camera parameters:
    [K, R, t] = read_xmp(filename);
    G = [R t; 0 0 0 1];
    %Full matrix
    ppm = K * [1 0 0 0
            0 1 0 0
            0 0 1 0] * G;
    %Project points to image
    [uv, vv] = proj(ppm,Points);
end