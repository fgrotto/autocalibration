close all;clear all;
addpath("./src");

directory = './dataset/zephyr_dante/';
load S.mat;

n_views = size(S);
for i = 1:5
    for j = 1:5
       if i == j 
           continue
       end
       I_left = [directory S{i,j}.name_view_i];
       I_right = [directory S{i,j}.name_view_j];
       
%        G = [S{i,j}.R_i S{i,j}.t_i; 0 0 0 1];
%        P1 = S{i,j}.K_i * [1 0 0 0
%                 0 1 0 0
%                 0 0 1 0] * G;
%             
%        G = [S{i,j}.R_j S{i,j}.t_j; 0 0 0 1];
%        P2 = S{i,j}.K_j * [1 0 0 0
%                 0 1 0 0
%                 0 0 1 0] * G;
%        
%        X = triang(P1,P2, [S{i,j}.uv_i(1:1) S{i,j}.vv_i(1:1)]', [S{i,j}.uv_j(1:1) S{i,j}.vv_j(1:1)]')
%        S{1,2}.points(1,:)

       if (size(S{i,j}.uv_i,1) > 50) && (size(S{i,j}.uv_j,1) > 50)
           leftP = [S{i,j}.uv_i(1:50) S{i,j}.vv_i(1:50)]';
           rightP = [S{i,j}.uv_j(1:50) S{i,j}.vv_j(1:50)]';

           F = fundamental_matrix(I_left, I_right, leftP, rightP);
           Fs(:,:,i,j) = F;
       end
    end
end

% Calculate initial estimation of the intrinsic camera parameters
image = imread([directory S{1,1}.name_view_i]);
info = imfinfo([directory S{1,1}.name_view_i]);
[H, W, C]= size(image);
SensorW=35;
Fmm=info.DigitalCamera.FocalLengthIn35mmFilm;
fp=(Fmm*W)/SensorW;
u_0=W/2;
v_0=H/2; 
K0=[fp 0 u_0; 0 fp v_0; 0 0 1];

% Calculate intrinsic camera parameters
K_vector = lsqnonlin(@(X) cost_medonca_cipolla(Fs, X),[K0(1,:) K0(2,2:3)],[],[], optimoptions('lsqnonlin','Display','off','Algorithm','levenberg-marquardt','TolX', 1e-10));

% intrinsic camera parameter in matrix form
K_cipolla = [K_vector(1) K_vector(2) K_vector(3); 
         0 K_vector(4) K_vector(5); 
         0 0 1];
     
K_auto = autocal(Fs,K0);

disp('Original');
disp(S{1,1}.K_i);
disp('Initial estimation');
disp(K0);
disp('Medonca cipolla: ');
disp(K_cipolla);
disp('Toolbox autocal: ');
disp(K_auto);

G = [S{1,2}.R_i S{1,2}.t_i; 0 0 0 1];
P1 = S{1,2}.K_i * [1 0 0 0
        0 1 0 0
        0 0 1 0] * G;

G = [S{1,2}.R_j S{1,2}.t_j; 0 0 0 1];
P2 = S{1,2}.K_j * [1 0 0 0
        0 1 0 0
        0 0 1 0] * G;

G21 = [K0\P2; 0 0 0 1]  * inv([ K0\P1; 0 0 0 1]);

for i = 1:5
    for j = 1:5
       if (size(S{i,j}.uv_i,1) > 50) && (size(S{i,j}.uv_j,1) > 50)
           leftP = [S{i,j}.uv_i(:) S{i,j}.vv_i(:)]';
           rightP = [S{i,j}.uv_j(:) S{i,j}.vv_j(:)]';

           % @TODO we are using the initial estimation use the estimated one
           % check the t21
           [R21,t21] = relative_lin(rightP, leftP, K0, K0);
           [R21,t21] = relative_nonlin(R21, t21, rightP, leftP, K0, K0);
           S{i,j}.R_ji = R21;
           S{i,j}.t_ji = t21;
       end
%       Just to check that we are getting the right R and t 
%       X = triang_lin_batch({P1, P2}, {leftP ,rightP});
%       [R1,t1] = exterior_lin(leftP,X,K0);
%       [R2,t2] = exterior_lin(rightP,X,K0);
     end
end

% Move every rotation to a rotation wrt of the first view
S{1,1}.R_new = eye(3);
S{1,1}.t_new = [0;0;0];
for i = 1:5
    for j = i:5
       if j ~= i
           S{i,j}.R_new = S{1, i}.R_ji * S{i,j}.R_ji;
           S{i,j}.t_new = S{1, i}.R_ji * S{i,j}.t_ji + S{1, i}.t_ji;
       else 
           S{i,j}.R_new = eye(3);
           S{i,j}.t_new = [0;0;0];
       end
    end
end


MyPoints = [];
for i = 1:5
    for j = (i+1):5
       if i == j 
           continue
       end
        G = [S{i,i}.R_new S{i,i}.t_new; 0 0 0 1];
        P1 = K0 * [1 0 0 0
                0 1 0 0
                0 0 1 0] * G;

        G = [S{i,j}.R_new S{i,j}.t_new; 0 0 0 1];
        P2 = K0 * [1 0 0 0
                0 1 0 0
                0 0 1 0] * G;
    
       if (size(S{i,j}.uv_i,1) > 50) && (size(S{i,j}.uv_j,1) > 50)
           leftP = [S{i,j}.uv_i(:) S{i,j}.vv_i(:)]';
           rightP = [S{i,j}.uv_j(:) S{i,j}.vv_j(:)]';
           X = triang_lin_batch({P1, P2}, {leftP ,rightP});
           MyPoints = [MyPoints; X'];
       end
    end
end

figure(1)
plot3(MyPoints(1:10:end, 1),MyPoints(1:10:end, 2), MyPoints(1:10:end,3), 'r.');
axis equal;
grid on;
