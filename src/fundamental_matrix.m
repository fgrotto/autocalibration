function F = fundamental_matrix(I_left, I_right, left_P, right_P)
  F0 = fund_lin(right_P,left_P);
  F_out = fund_nonlin(F0, right_P,left_P);
  
  F=F_out;
  
%   figure(1);
%   imshow(I_left);
%   figure(2);
%   imshow(I_right);
  
  % project points in the two images to see if they make sense
%   for i = 1:size(left_P,2)
%       figure(1); 
%       hold on;
%       plot(left_P(1,i),left_P(2,i),'r*');
% 
%       figure(2);
%       hold on;
%       plot(right_P(1,i),right_P(2,i),'r*');
%       disp('Done');
%   end
  
%   %Draw epipolar lines:
%   disp('Draw 3 epipolar lines: ');
%   [m n c] = size(I_left);
%   list =['c' 'b' 'g'];
%   
%   for i=1:3
%         % Clicking a point on the left image:
%         figure(1);    
%         [left_x, left_y] = ginput(1);
%         hold on;
%         plot(left_x,left_y,'r*');
% 
%         % Finding the epipolar line on the right image:
%         left_P = [left_x; left_y; 1];
% 
%         right_P = F*left_P;
% 
%         right_epipolar_x=1:n;
%         % Using the eqn of line: ax+by+c=0; y = (-c-ax)/b
%         right_epipolar_y=(-right_P(3)-right_P(1)*right_epipolar_x)/right_P(2);
%         figure(2);
%         hold on;
%         plot(right_epipolar_x,right_epipolar_y,list(mod(i,8)));
%   end
end