function G = absolute(X,Y)

% ABSOLUTE solve absolute orientation 
% i.e. compute rigid transformation between two 3D point sets X and Y,
% When G is applied to Y, it brings Y on X. Entries containing NaN in X
% are discarded together with the corresponding entries in Y (used for 
% discarding correspondences)
% Algorithm ref.: Kanatani
%
% Author: A. Fusiello, 1998

% discard NaN entries in X and correspondingly in Y
i = find(~isnan(X));
X = reshape(X(i),length(X(i))/3,3);
Y = reshape(Y(i),length(Y(i))/3,3);

dime = size(Y,1);

% compute centroids
cm = sum(Y,1)./dime;
cd = sum(X,1)./dime;

% subtract centroids
Yb = rigid([0,0,0, -cm]',Y);
Xb  = rigid([0,0,0, -cd]',X);

% compute rotation
K =  Xb' * Yb;
[U,D,V]=svd(K);
S = diag([1,1,det(U*V')]);
R = U*S*V';

%  compute traslation
t = cd' - R * cm';

% rigid transformation
G =  [R, t
    0 0 0 1];
