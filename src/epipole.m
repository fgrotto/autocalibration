function e = epipole(F)
    % Extract epipole from fundamental matrix via SVD
    [~, ~, V]  = svd(F');
    x = V(:,end);
    e = [ 0 -x(3) x(2);
         x(3) 0 -x(1); 
        -x(2) x(1) 0];
end