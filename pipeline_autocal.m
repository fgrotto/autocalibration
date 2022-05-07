close all;clear all;
addpath("./src");

directory = './dataset/zephyr_dante/';
load S.mat;

% Extract dimension of dataset
n_views = size(S,1);
% Consider only this number of points per image
num_points = 30;
% Consider only a limited number of pair of images
image_indexes = randi(n_views,20,1); % TODO maybe 1:1:20

% Compute fundantal matrixes considering the provided indexes and
% the number of points for each pair of images
Fs = compute_f(S, directory, image_indexes, num_points);

% Calculate initial estimation of the intrinsic camera parameters
K0 = compute_k0(S, directory);

% Compute intrinsic parameters with medonca cipolla
K_cipolla = compute_mc(Fs,K0);    
K_cipolla_toolbox = autocal(Fs,K0);

disp('Estimation of intrinsic parameters');
disp('Original');
disp(S{1,1}.K_i);
disp('Initial estimation');
disp(K0);
disp('Medonca cipolla: ');
disp(K_cipolla);
disp('Medonca cipolla (Toolbox): ');
disp(K_cipolla_toolbox);


%Project original points
Points = [];
for i = 1:1
    for j = (i+1):2
        Points = [Points; S{i,j}.points];
    end
end


% Plot points with the mesh
% figure(3)
% [Tri,Pts] = plyread([directory 'Mesh.ply'],'tri');
% trisurf(Tri,Pts(:,1),Pts(:,2),Pts(:,3)); 
% colormap(gray); axis equal;
% hold on
% plot3(Points(1:10:end, 1),Points(1:10:end, 2), Points(1:10:end,3), 'r.');


for i = 1:5
    for j = (i+1):5
       if (size(S{i,j}.uv_i,1) > 50) && (size(S{i,j}.uv_j,1) > 50)
           leftP = [S{i,j}.uv_i(1:50) S{i,j}.vv_i(1:50)]';
           rightP = [S{i,j}.uv_j(1:50) S{i,j}.vv_j(1:50)]';
           [R21,t21] = relative_lin(rightP, leftP, K_cipolla, K_cipolla);
           [R21,t21] = relative_nonlin(R21, t21, rightP, leftP, K_cipolla, K_cipolla);
           S{i,j}.R_ji = R21;
           S{i,j}.t_ji = t21;
       end
    end
end

MyPoints = [];
S{1,1}.R_ji = eye(3);
S{1,1}.t_ji = [0;0;0];
S{1,1}.R_new = eye(3);
S{1,1}.t_new = [0;0;0];

% Problem in case of missing steps we can't get the correct points from
% view 1, esample missing 1-8 will prevent to build 8-i views
for i = 1:1
    for j = (i+1):3
        S{i,j}.R_new = S{1, i}.R_ji * S{i,j}.R_ji;
        S{i,j}.t_new = S{1, i}.R_ji * S{i,j}.t_ji + S{1, i}.t_ji;
        
        G = [S{1,i}.R_new S{1,i}.t_new; 0 0 0 1];
        P1 = K_cipolla * [1 0 0 0; 0 1 0 0; 0 0 1 0] * G;

        G = [S{i,j}.R_new S{i,j}.t_new; 0 0 0 1];
        P2 = K_cipolla * [1 0 0 0; 0 1 0 0; 0 0 1 0] * G;

       if (size(S{i,j}.uv_i,1) > 50) && (size(S{i,j}.uv_j,1) > 50)
           leftP = [S{i,j}.uv_i(:) S{i,j}.vv_i(:)]';
           rightP = [S{i,j}.uv_j(:) S{i,j}.vv_j(:)]';
           X = triang_lin_batch({P1, P2}, {leftP ,rightP});
           MyPoints = [MyPoints; X'];
       end
    end
end

% Get two clouds of points
cloud_points1 = Points;
cloud_points2 = MyPoints;

% Apply ICP algorithm and calculate the rigid transformation
model = cloud_points1;
data = cloud_points2;
Gi = icp(model,data);
datafinal = rigid(Gi,data);

% Display the original mesh, the original 3d points and the estimated
% points with our procedure
figure(3)
[Tri,Pts] = plyread([directory 'Mesh.ply'],'tri');
trisurf(Tri,Pts(:,1),Pts(:,2),Pts(:,3)); 
colormap(gray); axis equal;
hold on
plot3(model(:,1), model(:,2), model(:,3), '.b');
hold on
plot3(datafinal(:,1), datafinal(:,2), datafinal(:,3), '.r');
hold on
grid on
axis equal
