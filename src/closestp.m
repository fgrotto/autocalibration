function  [res,modelcp] = closestp(data,model)

% CLOSESTP     compute closest points
% modelcp(i,:) are the coordinates of the closest point in the model set
% of the i-th point in the data set
% res is the residual, i.e., the average distance between closest points
% 
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

res = mean(mindist);
