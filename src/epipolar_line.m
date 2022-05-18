function lines = epipolar_line(F,pts)
%   The output is a M-by-3 matrix where each row has the format of [A,B,C]
%   which defines a line as A * x + B * y + C = 0. M is the number of
%   lines.

    lines = [pts, ones(size(pts, 1), 1)] * F';
end
