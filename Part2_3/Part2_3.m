% Extended Kalman Filter for tricycle Mobile Robot type 2
clear all
close all

load measdat6
load contrlsig
 
G = 6;  % Group 6
dt = 0.05;  % Time step;
runTime = 40; % Total run time
step = floor(runTime/dt); %Step
v = 0.1; % Velocity
alpha = 0.2; % Steerinng angle
d = sqrt(2.5+0.1*G); %Axial distance

% Landmarks
x_LM1 = 3;
y_LM1 = 4;

x_LM2 = 4;
y_LM2 = 4;

% Initial error covarience matrix
P(:,:,1)=[0 0 0;0 0 0;0 0 0];     
                                            
% Initial location
x(1) = 0;  
y(1) = 0;
theta(1) = 0;

% Initial actual location  
xa(1) = 0;  
ya(1) = 0;
thetaa(1) = 0;

% Initial prediction location 
xp(1) = 0;  
yp(1) = 0;
thetap(1) = 0;

% Variance of the input noise.
Q=[(2.5*10^-3) 0; 0 (3.6*10^-3)];    

% Variance of the measurement noise.
R=[(10^-6) 0 0 0;...
    0 (10^-6) 0 0;...
    0 0 (7.62*10^-5) 0;...
    0 0 0 (7.62*10^-5)];  

%Input control
v_a = inp(1,:);
alpha_a = inp(2,:);

z(:,:) = me(:,:);

for k=2:step
    % Actual value
    xa(k) = xa(k-1) + v_a(k-1)*cos(thetaa(k-1))*dt;
    ya(k) = ya(k-1) + v_a(k-1)*sin(thetaa(k-1))*dt;
    thetaa(k) = thetaa(k-1) + v_a(k-1)*tan(alpha_a(k-1))/d*dt;

    % Prediction based on system model only (no Kalman Filter)
    xp(k) = xp(k-1) + v*cos(thetap(k-1))*dt;
    yp(k) = yp(k-1) + v*sin(thetap(k-1))*dt;
    thetap(k) = thetap(k-1) + (v*tan(alpha))/d*dt;
   
    % Kalman filter
    A = [1 0 -v*sin(theta(k-1))*dt;...
         0 1 v*cos(theta(k-1))*dt;...
         0 0 1];
    
    W = [cos(theta(k-1))*dt 0;...
         sin(theta(k-1))*dt 0;...
         tan(alpha)*dt/d v*dt/(cos(alpha)*cos(alpha)*d)];
     
    Pe(:,:,k)=A*P(:,:,k-1)*A' + W*Q*W';
     
    H = [(x(k-1)-x_LM1)/sqrt((x(k-1)-x_LM1)^2 + (y(k-1)-y_LM1)^2)...
       (y(k-1)-y_LM1)/sqrt((x(k-1)-x_LM1)^2 + (y(k-1)-y_LM1)^2) 0;...
       (x(k-1)-x_LM2)/sqrt((x(k-1)-x_LM2)^2 + (y(k-1)-y_LM2)^2)...
       (y(k-1)-y_LM2)/sqrt((x(k-1)-x_LM2)^2 + (y(k-1)-y_LM2)^2) 0;...
       (y_LM1-y(k-1))/((x(k-1)-x_LM1)^2 + (y(k-1)-y_LM1)^2)...
       (x_LM1-x(k-1))/((x(k-1)-x_LM1)^2 + (y(k-1)-y_LM1)^2) -1;...
       (y_LM2-y(k-1))/((x(k-1)-x_LM1)^2 + (y(k-1)-y_LM1)^2)...
       (x_LM2-x(k-1))/((x(k-1)-x_LM1)^2 + (y(k-1)-y_LM1)^2) -1];
    
   % Kalman gain
    K(:,:,k)=Pe(:,:,k)*H'/(H*Pe(:,:,k)*H' + R);  

    
    X_pre(:,k) = [xa(k-1) + v*cos(thetaa(k-1))*dt;
                  ya(k-1) + v*sin(thetaa(k-1))*dt;
                  thetaa(k-1) + v*tan(alpha)/d*dt];
    
    rA = sqrt((X_pre(1,k-1)-x_LM1).*(X_pre(1,k-1)-x_LM1) + (X_pre(2,k-1)-y_LM1).*(X_pre(2,k-1)-y_LM1));
    rB = sqrt((X_pre(1,k-1)-x_LM2).*(X_pre(1,k-1)-x_LM2) + (X_pre(2,k-1)-y_LM2).*(X_pre(2,k-1)-y_LM2));
    phi_A = atan2((y_LM1-X_pre(2,k-1)),(x_LM1-X_pre(1,k-1))) - X_pre(3,k-1);
    phi_B = atan2((y_LM2-X_pre(2,k-1)),(x_LM2-X_pre(1,k-1))) - X_pre(3,k-1);
    z_pre(:,k) = [rA; rB; phi_A; phi_B];

    % Kalman Estimated location
    X = X_pre(:,k) + K(:,:,k)*(z(:,k)- z_pre(:,k));  
    x(k) = X(1);
    y(k) = X(2);
    theta(k) = X(3);
    
    % Estimation covariance error
    P(:,:,k)=(eye(3)-K(:,:,k)*H)*Pe(:,:,k); 
end

% Plot the results
time=0:step-1;
figure(1)
plot(x,y,xp,yp,xa,ya,'--')
xlabel('X (m)');
ylabel('Y (m)');
legend('Kalman','Prediction','Actual','Location','SouthEast' )
title('Actual and Estimated Trajectories')
grid on;

% Creating the zoom 
ax = axes;
set(ax,'units','normalized','position',[0.2,0.6,0.3,0.3])
box(ax,'on')
hold on
plot(x,y,xp,yp,xa,ya,'--')
set(ax,'xlim',[3.3,3.4],'ylim',[0.65,0.69])
grid on;

figure(2)
plot(time,x-xa,time,xp-xa,'--')
legend('Kalman error','Prediction error')
xlabel('Time (step)');
ylabel('Error (m)');
title('Error in X')
grid

figure(3)
plot(time,y-ya,time,yp-ya,'--')
xlabel('Time (step)');
ylabel('Error (m)');
legend('Kalman error','Prediction error')
title('Error in Y')
grid

figure(4)
plot(time,theta-thetaa,time,thetap-thetaa,'--')
xlabel('Time (step)');
ylabel('Error (rad)');
legend('Kalman error','Prediction error')
title('Error in Theta')
grid