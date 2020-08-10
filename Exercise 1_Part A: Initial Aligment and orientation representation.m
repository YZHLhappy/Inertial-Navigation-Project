%2020.04.28 Yang
clear all
load('imubsp.mat')
format long g

%Find initial orientation for attitude computation
s = size(imudata);
z = s(1);
%import the data (timestamp,fx,fy,fz,wx,wy,wz)
timestamp = [];
f_b_x = [];
f_b_y = [];
f_b_z = [];
w_b_x = [];
w_b_y = [];
w_b_z = [];

for k = 1:z
    timestamp(k,:) = imudata(k,1);
    f_b_x(k,:) = imudata(k,5);
    f_b_y(k,:) = imudata(k,6);
    f_b_z(k,:) = imudata(k,7);
    w_b_x(k,:) = imudata(k,2);
    w_b_y(k,:) = imudata(k,3);
    w_b_z(k,:) = imudata(k,4);
end

f_b_x_mean = mean(f_b_x(:));
f_b_y_mean = mean(f_b_y(:));
f_b_z_mean = mean(f_b_z(:));

w_b_x_mean = mean(w_b_x(:));
w_b_y_mean = mean(w_b_y(:));
w_b_z_mean = mean(w_b_z(:));
w_xyz = [w_b_x_mean;w_b_y_mean;w_b_z_mean];

fai0 = atan2(-f_b_y_mean,-f_b_z_mean);
theta0 = atan(-f_b_x_mean/sqrt(f_b_y_mean^2+f_b_z_mean^2));

c2 = [cos(theta0) 0 sin(theta0); 0 1 0; -sin(theta0) 0 cos(theta0)];
c1 = [1 0 0; 0 cos(fai0) -sin(fai0); 0 sin(fai0) cos(fai0)];

w_b_n_xyz = c1 * c2 * w_xyz; 
w_b_n_x = w_b_n_xyz(1);
w_b_n_y = w_b_n_xyz(2);

psi0 = atan(w_b_n_y/w_b_n_x);

c3 = [cos(psi0) -sin(psi0) 0;sin(psi0) cos(psi0) 0;0 0 1];

%fai0,theta0,psi0
fai0_deg=rad2deg(fai0);
theta0_deg=rad2deg(theta0);
psi0_deg=rad2deg(psi0);

%c_n->b
c_n_to_b = c1*c2*c3;

%g
acc = [f_b_x_mean;f_b_y_mean;f_b_z_mean];
g = sqrt(acc(1)^2+acc(2)^2+acc(3)^2);
%w_e
b = sqrt(w_xyz(1)^2+w_xyz(2)^2+w_xyz(3)^2);
b_deg_h = rad2deg(b)*3600;

%trace of the matrix,u,delta
d = trace(c_n_to_b);
delta = acos((d-1)/2);
vector_delta = [(c_n_to_b(3,2)-c_n_to_b(2,3));
                (c_n_to_b(1,3)-c_n_to_b(3,1));
                (c_n_to_b(2,1)-c_n_to_b(1,2))];
u = vector_delta/norm(vector_delta);            

%quaternion,a,b,c,d
Q_a = (1/2)*sqrt(1+trace(c_n_to_b));
Q_part = (1/(4*Q_a))*vector_delta;
Q_b = Q_part(1);
Q_c = Q_part(2);
Q_d = Q_part(3);
Q = [Q_a;Q_b;Q_c;Q_d];
