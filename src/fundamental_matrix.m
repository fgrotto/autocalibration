function [F, left_in, right_in, inliers] = fundamental_matrix(I_left, I_right, left_P, right_P, num_points)
    %   image_indexes = randi(size(left_P),num_points,1);
    %   left_P = left_P(:, image_indexes);
    %   right_P = right_P(:, image_indexes);

    %   F0 = fund_lin(right_P,left_P);
    %   F_out = fund_nonlin(F0, right_P,left_P);

    [F_out, inliers] = estimateFundamentalMatrix(left_P', right_P', 'Method', 'RANSAC', 'NumTrials', 1000, 'DistanceThreshold', 1e-4);
    F = F_out;

    matchedPoints1 = left_P';
    matchedPoints2 = right_P';

    left_in = matchedPoints1(inliers, :)';
    right_in = matchedPoints2(inliers, :)';
%     I1 = imread(I_left);
%     figure;
%     subplot(211);
%     imshow(I1);
%     title('Inliers and Epipolar Lines in First Image'); hold on;
%     plot(matchedPoints1(inliers, 1), matchedPoints1(inliers, 2), 'go')
%     epiLines = epipolarLine(F', matchedPoints2(inliers, :));
%     points = lineToBorderPoints(epiLines, size(I1));
%     line(points(:, [1, 3])', points(:, [2, 4])');
%     I2 = imread(I_right);
%     subplot(212);
%     imshow(I2);
%     title('Inliers and Epipolar Lines in Second Image'); hold on;
%     plot(matchedPoints2(inliers, 1), matchedPoints2(inliers, 2), 'go')
%     epiLines = epipolarLine(F, matchedPoints1(inliers, :));
%     points = lineToBorderPoints(epiLines, size(I2));
%     line(points(:, [1, 3])', points(:, [2, 4])');
%     truesize;
end
