function  [res,modelcp] = closestp_x84(data,model)

% CLOSESTP     compute closest points
% modelcp(i,:) are the coordinates of the closest point in the model set
% of the i-th point in the data set
% res is the residual, i.e., the average distance between closest points
% The X84 strategy is applied to discard outliers
%
% See also ICP

modelcp = ones(size(data));
mindist = inf*ones(size(data,1),1);

for i = 1:size(data,1)
    for j = 1:size(model,1)    
        d = norm(model(j,:) - data(i,:));
        if d < mindist(i)
            mindist(i)=d;
            modelcp(i,:) =  model(j,:);
        end
    end
end

% apply X84
location = median(mindist);
scale = 5.2 * median(abs(mindist-location));

% set points to NaN in order to discard them
I = find(abs(mindist-location) > scale);
modelcp(I,:)=NaN;

% compute average distance for inliers
J = find(abs(mindist-location) <= scale);
res = mean(mindist(J));
