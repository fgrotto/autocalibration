function D  = rigid(G,M)

% RIGID    apply 3D rigid transformation G to a point set M
% the transformation can be a 4x4 homogeneous matrix or a 
% vector whose first 3 components are the rotation angles
% around the x,y,z axis resp. and the last 3 are the translation
% vector.
%
% See also EUL

if size(G,2) == 1
   % is a vector
   G = [eul(G(1:3)) G(4:6)];
else
   % is a matrix
   G = G(1:3,:);
end

HM =[M ones(size(M,1),1)]';
D = (G*HM)';
