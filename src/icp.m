function G = icp(model, data, maxiter)

% ICP   Iterative Closet Point algorithm.
% Align a cloud of points (data) to another (model) 
% Return a 4x4 homogeneous matrix representing a rigid trasform
% When G is applied to data, it brings data onto model.
% 
% See also RIGID

% Author: A. Fusiello, 2003

if (nargin < 3)
   maxiter = 200;
end

%eps = 0.0001;
eps = 0.00000001;

G=eye(4);
res = inf;
prevres = 0;
i=0; % iterations counter
while ((abs(res-prevres)> eps) && i < maxiter)
    i = i+1;
    prevres = res;
    % apply current transformation bringing data onto model 
    dataREG = rigid(G,data);
    [mindist,modelCP] = closestp(dataREG,model);
 
    res = mean(mindist);

    % compute incremental tranformation
    GI = absolute(modelCP,dataREG);
    G = GI * G;
end

fprintf('Iterations: %d\n', i);
