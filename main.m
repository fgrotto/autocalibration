clear all;
close all;
warning('off','all')

addpath("./src");
load('dataset/autocalib.mat');

%% Autocalibration from Computer Vision Toolbox and medonca-cipolla method
% initial K approximation
K_approx = A;
disp('Initial K approximation: ')
disp(K_approx)

% Compute medonca-cipolla method
K_vector = lsqnonlin(@(X) cost_medonca_cipolla(Fs, X),[K_approx(1,:) K_approx(2,2:3)],[],[], optimoptions('lsqnonlin','Display','off','Algorithm','levenberg-marquardt','TolX', 1e-10));

% intrinsic camera parameter in matrix form
K = [K_vector(1) K_vector(2) K_vector(3); 
         0 K_vector(4) K_vector(5); 
         0 0 1];
     
disp('Estimated internal camera parameters (Medonca-Cipolla): ');
disp(K);

disp('Estimated internal camera parameters (Autocalibration Toolbox): ');
disp(autocal(Fs,K_approx));

%% Practical Autocalibration (Riccardo Gherardi and Andrea Fusiello)

[best_H, best_focal_length, best_cost] = practical_autocal(1,10,PPM);

disp('Best H');
disp(best_H);

disp('Best focal length');
disp(best_focal_length);

disp('Best cost');
disp(best_cost);