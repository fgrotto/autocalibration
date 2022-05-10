close all;clear all;
addpath("./src");

directory = './dataset/zephyr_dante/';
load S.mat;

normF = @(x) norm(x,'fro');
disp('Pipeline for autocalibration starting from S datastructure');

% Extract dimension of dataset
n_views = size(S,1);
% Consider only this number of points per image
num_points = 30;
% Consider only a limited number of pair of images
image_indexes = randi(n_views,20,1); % TODO maybe 1:1:20

% Compute fundantal matrixes considering the provided indexes and
% the number of points for each pair of images
disp('Computing fundamental matrices');
Fs = compute_f(S, directory, image_indexes, num_points);

% Calculate initial estimation of the intrinsic camera parameters
K0 = compute_k0(S, directory);

% Compute intrinsic parameters with medonca cipolla and kruppas method
K_cipolla = compute_mc(Fs,K0);   
K_kruppas = compute_kruppas(Fs,K0);
K_cipolla_toolbox = autocal(Fs,K0);

disp('Estimating intrinsic parameters');
disp('Original');
disp(S{1,1}.K_i);
disp('Initial estimation');
disp(K0);
disp('Medonca cipolla: ');
disp(K_cipolla);
disp('Kruppas method: ');
disp(K_kruppas);
disp('Medonca cipolla (Toolbox): ');
disp(K_cipolla_toolbox);


% Project original points and collect them in a proper vector
Points = [];
for i = 1:1
    for j = (i+1):3
        if (size(S{i,j}.uv_i,1) > 50) && (size(S{i,j}.uv_j,1) > 50)
            Points = [Points; S{i,j}.points];
        end
    end
end

disp('Computing relative orientations of views');
for i = 1:5
    for j = (i+1):5
       if (size(S{i,j}.uv_i,1) > 50) && (size(S{i,j}.uv_j,1) > 50)
           % Get around 50 points per view
           leftP = [S{i,j}.uv_i(1:50) S{i,j}.vv_i(1:50)]';
           rightP = [S{i,j}.uv_j(1:50) S{i,j}.vv_j(1:50)]';
           
           % Compute original transformation matrix to get the actual error
           G1 = [S{i,j}.R_i S{1,i}.t_i; 0 0 0 1];
           G2 = [S{i,j}.R_j S{i,j}.t_j; 0 0 0 1];
           G21 = G2*inv(G1);
           
           [R21,t21] = relative_lin(rightP, leftP, K_cipolla, K_cipolla);
           fprintf('Relative linear SO3 views (%0.2g, %0.2g) error:\t %0.5g \n', i,j, normF(R21 - G21(1:3,1:3)));
           [R21,t21] = relative_nonlin(R21, t21, rightP, leftP, K_cipolla, K_cipolla);
           fprintf('Relative nonlin SO3 views (%0.2g, %0.2g) error:\t %0.5g \n', i,j, normF(R21 - G21(1:3,1:3) ));
           if normF(R21 - G21(1:3,1:3)) > 1
               continue
           end
           % Save the obtained rotation and translation of j wrt of i
           S{i,j}.R_ji = R21;
           S{i,j}.t_ji = t21;
           
           % Triangulate the points with the obtained rotation and translation 
           x1 = leftP;
           x2 = rightP;
           X_model = triang_lin_batch({K_cipolla*[eye(3),zeros(3,1)], K_cipolla*[R21,t21]}, {x1,x2});
           X = S{i,j}.points(1:50,:)';
           
           % Apply opa to the obtained 3D points and the original ones to
           % recover the actual rotation, translation and scale
           [R,t,s] = opa(X(:,1:6),X_model(:,1:6));
           X_obj = s*(R*X_model + t*ones(1,size(X,2)));
           fprintf('Relative triang error views (%0.2g, %0.2g):\t %0.5g \n',i,j, rmse(X(:)-X_obj(:)));

           figure;
           plot3(X(1,:), X(2,:), X(3,:), 'or'); hold on
           plot3(X_obj(1,:), X_obj(2,:), X_obj(3,:),'+b');
           title(strcat('Relative orientation views (',num2str(i),',',num2str(j),')'));
       end
    end
end

% Concatenate all the points
disp('Computing final points by view concatenation');
MyPoints = [];
S{1,1}.R_ji = eye(3);
S{1,1}.t_ji = [0;0;0];
S{1,1}.R_new = eye(3);
S{1,1}.t_new = [0;0;0];

% Problem in case of missing steps we can't get the correct points from
% view 1, example missing 1-8 will prevent to build 8-i views
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
           
           % Triangulate the points with the obtained rotation and translation
           x1 = leftP;
           x2 = rightP;
           X_model = triang_lin_batch({P1, P2}, {leftP ,rightP});
           X = S{i,j}.points';
           
           % Apply opa to the obtained 3D points and the original ones to
           % recover the actual rotation, translation and scale
           [R,t,s] = opa(X(:,1:6),X_model(:,1:6));
           X_obj = s*(R*X_model + t*ones(1,size(X,2)));
           fprintf('Relative triang error views (%0.2g, %0.2g):\t\t %0.5g \n', i,j, rmse(X(:)-X_obj(:)));
            
           % Plot figure to show the obtained points
           figure;
           plot3(X(1,:), X(2,:), X(3,:), 'or'); hold on
           plot3(X_obj(1,:), X_obj(2,:), X_obj(3,:),'+b');
           title(strcat('Relative orientation views (',num2str(i),',',num2str(j),')'));
           
           MyPoints = [MyPoints X_obj];
       end
    end
end

disp('Computing icp with the original model and the obtained points');
% Get two clouds of points
cloud_points1 = Points;
cloud_points2 = MyPoints';

% Apply ICP algorithm and calculate the rigid transformation
model = cloud_points1;
data = cloud_points2;
Gi = icp(model,data);
datafinal = rigid(Gi,data);
disp('Pipeline completed');

% Display the original mesh, the original 3d points and the estimated
% points with our procedure
figure;
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
