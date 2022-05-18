% Simple script to obtain a datastructure for automatic calibration
% procedure. The idea is to create a datastructure for each pair of images
% containing the 3D point and the two corrisponding points on image i,j (we
% will project only points that as visible in both images

close all;clear all;
addpath("./src");
directory = './dataset/zephyr_dante/';

% Read and show sparse point clouds: (there are many noisy points, it is good to subsample the cloud) 
[Points] = plyread([directory 'SamPointCloud.ply']);
X=[Points.vertex.x Points.vertex.y Points.vertex.z];
figure(1)
plot3(X(1:10:end, 1),X(1:10:end, 2), X(1:10:end,3), 'r.');
axis equal;
grid on;

% Read and show mesh for points
figure(2)
[Tri,Pts] = plyread([directory 'Mesh.ply'],'tri');
trisurf(Tri,Pts(:,1),Pts(:,2),Pts(:,3)); 
colormap(gray); axis equal;
hold on
plot3(X(1:10:end, 1),X(1:10:end, 2), X(1:10:end,3), 'r.');
 
% Open and display the title of the visibility points
fid=fopen([directory 'Visibility.txt'], 'r');
line_ex = fgetl(fid);
disp(line_ex);

nviews= str2num(fgetl(fid));
disp(['Number of cameras\views: ' num2str(nviews)])

% Load all images and cameras
for i=1:nviews
    string=fgetl(fid);
    name_view_i=string(end-11:end);
    
    % Count the number of visible points and initialize the 
    % array with the values read from the visibility file
    npoint=str2num(fgetl(fid));
    VisPoints=zeros(npoint,3);
    for p=1:npoint
        VisPoints (p,:)=str2num(fgetl(fid));
        VisPoints(p, 1)= VisPoints(p, 1)+1;
    end
    
    % Extract visible points from image i
    Xvis_i=X(VisPoints(:,1),:);
    
    tmp{i}.points = Xvis_i;
    tmp{i}.name_view = name_view_i;
end

for i = 1:size(tmp,2)
    for j = 1:size(tmp,2)
        Xvis_i = tmp{i}.points;
        Xvis_j = tmp{j}.points;
        name_view_i = tmp{i}.name_view;
        name_view_j = tmp{j}.name_view;
        
        % Extract intersection from i,j 3d points defined by row
        [CommonPointsOriginal,ia,ib] = intersect(Xvis_i,Xvis_j,'rows');
        
        % Extract points close to the centroid
        CommonPoints = extractPointsUsingCentroid(CommonPointsOriginal, 3);
        if isempty(CommonPoints)
           CommonPoints = CommonPointsOriginal;
        end
        
        % Get points projections to image i
        [uv_i, vv_i, K_i, R_i, t_i] = points_to_image([directory name_view_i(1:end-3) 'xmp'], CommonPoints); 
        
        % Get points projections to image j
        [uv_j, vv_j, K_j, R_j, t_j] = points_to_image([directory name_view_j(1:end-3) 'xmp'], CommonPoints);
        
        S{i,j}.points = CommonPoints;
        S{i,j}.uv_i = uv_i;
        S{i,j}.vv_i = vv_i;
        S{i,j}.uv_j = uv_j;
        S{i,j}.vv_j = vv_j;
        S{i,j}.name_view_i = name_view_i;
        S{i,j}.name_view_j = name_view_j;
        S{i,j}.K_i = K_i;
        S{i,j}.R_i = R_i;
        S{i,j}.t_i = t_i;
        S{i,j}.K_j = K_j;
        S{i,j}.R_j = R_j;
        S{i,j}.t_j = t_j;
    end
end

% Save the structure for the pipeline
save('S.mat', 'S');

