function a  = ieul(R)

% IEUL: compute euler angles (yaw,pitch,roll) givena rotation matrix 

phi = atan2(R(2,1),R(1,1)); 
theta = asin(- R(3,1));
psi =  atan2(R(3,2),R(3,3));

a(1)=psi;
a(2)=theta;
a(3)=phi;
