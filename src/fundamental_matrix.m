function [F] = fundamental_matrix(I_left, I_right, left_P, right_P)
    F0 = fund_lin(right_P,left_P);
    F  = fund_nonlin(F0, right_P,left_P);
    
%     % Select only 20 points for epipolar lines
%     matchedPoints1 = left_P(:,30:50)';
%     matchedPoints2 = right_P(:,30:50)';
%     
%     I1 = imread(I_left);
%     figure;
%     subplot(211);
%     imshow(I1);
%     title('Epipolar Lines in First Image'); hold on;
%     plot(matchedPoints1(:, 1), matchedPoints1(:, 2), 'go')
%     epiLines = epipolarLine(F', matchedPoints2(:, :));
%     points = lineToBorderPoints(epiLines, size(I1));
%     line(points(:, [1, 3])', points(:, [2, 4])');
%     I2 = imread(I_right);
%     subplot(212);
%     imshow(I2);
%     title('Epipolar Lines in Second Image'); hold on;
%     plot(matchedPoints2(:, 1), matchedPoints2(:, 2), 'go')
%     epiLines = epipolarLine(F, matchedPoints1(:, :));
%     points = lineToBorderPoints(epiLines, size(I2));
%     line(points(:, [1, 3])', points(:, [2, 4])');
%     truesize;
end
