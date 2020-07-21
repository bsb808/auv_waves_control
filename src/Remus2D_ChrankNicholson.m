%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   2D Dynamic analysis of REMUS in xz-plane
%
%         17July2020
%
%  This program comuputes the dynamic motion of the REMUS underwater vehicles
%   subjected to various hydrodynamic loads in the xz-plane.
%    
%   [Amat]{acc}+[Bmat]{ve}={extLoad} is soleved using the
%   Crank-Nicholson technique which has the tuncation order (dt)^2
%   Here, {vel}={u; w; q} column vector
%         {acc}=(udot; wdot; qdot} column vector
%         {extLoad} is the force and momemnt indpendent of u, w, and q.
%
%------------------------------------------------------------------------
%  State varaibles
%   u;  linear velocity along the x-axis
%   w:  linear velocity along the z-axis
%   q:  angular velocity about the y-axis
%
%   udot:   linear acceleration along the x-axis
%   vdot:   linear acceleration along the z-axis
%   qdot:   angular avveleration about the y-axis
%
%   x_position: position along the x-axis
%   z_position: position along the z-axis
%   q_position: angular position about y-axis
%
%   Fx: linear force along the x-axis
%   Fy: linear force along the z-axis
%   My: moment about the y-axis
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc

m=30;       % mass in kg
Iyy=3.45;   % mass moment of interia about y-axis in kg-m^2

xg=0;       % center of gravisty wrt center of bouyance in x-axis 
zg=1.96e-2; % center of gravisty wrt center of bouyance in z-axis (m)

Xuu=-1.62;          % Axial flow drag (Kg/m)
Xudot=-9.30e-001;	% Added Mass (kg)
Xwq=-3.55e+001;     % Added Mass Cross-term (kg/rad)
Xqq=1.93e+000;      % Added Mass Cross-term (kg-m/rad)
Xprop=+3.86e+000;   % Propeller Thrust (N)
Zww=-1.31e+002;     % Cross-flow Drag (kg/m)
Zqq=-6.32e-001;     % Cross-flow Drag (kg· m/rad^2)
Zuw=-2.86e+001;  	% Body Lift Force and Fin Lift (Kg/m)
Zwdot=-3.55e+001;   % Added Mass (kg)
Zqdot=-1.93e+000;   % Added Mass (kg• m/rad)
Zuq=-5.22e+000;     % Added Mass Cross-term and Fin Lift (kg/m)
Zuud=-9.64e+000;    % Fin Lift Force (kg/(m · rad))
Mww=+3.18e+000;     % Crossflow Drag (kg)
Mqq=-9.40e+000;     % Crossflow Drag (kg· m^2 / rad^2)
Muw=+2.40e+001;     % Body and Fin Lift and Munk Moment (kg)
Mwdot=-1.93e+000;   % Added mass (kg-m)
Mqdot=-4.88e+000;   % Added mass (kg-m^2 /rad)
Muq=-2.00e+000;     % Added Mass Cross Term and Fin Lift (kg-m/rad)
Muud=-6.15e+000;    % Fin Lift Moment (kg/rad)

del=0.0;           % fin angle

Amat=zeros(3,3);    % initialization of matrice Amat
Bmat=zeros(3,3);    % initialization of matrice Bmat
exLoad=zeros(3,1);  % initialization of column vector exLoad
Cmat=zeros(3,3);    % temporary matrix
FF=zeros(3,3);      % temporary vector

dt=0.001;      % time step size (sec)
ntime=2000;      % no of time steps
nit=5000;     % no of iteration for convergence
tol=0.001;   % tolerance for convergence

u=zeros(ntime+1,1);     % velocity in x-axis
w=zeros(ntime+1,1);     % velocity in z-axis
q=zeros(ntime+1,1);     % angular velocity 
x_position=zeros(ntime+1,1);        % position in x-axis
z_position=zeros(ntime+1,1);        % position in z-axis
q_position=zeros(ntime+1,1);        % angular position
time=zeros(ntime+1,1);     % time

NoIt=zeros(ntime+1,1);       % array to store no. of iteration for convergency
ConvError=zeros(ntime+1,1);       % array to store the error ratio for convergency

u(1)=0;     % inital x-velocity
v(1)=0;     % inital z-velocity
q(1)=0;     % inital angular velocity about y-axis

	
for itime=1:ntime   % time increment loop
    
    time(itime+1)=time(itime)+dt;
    
    up=u(itime);     % current velocity in x-axis
    wp=w(itime);     % current velocity in z-axis
    qp=q(itime);     % current angular velocity about z-axis

    it=0;
    ratio=1;     % initial value for covenrgence

    while (it <= nit && ratio > tol)
    it=it+1;    
    solo=[up; wp; qp];      % Assumed solution
    
    Amat(1,1)=m-Xudot;  % Mtarix A
    Amat(1,3)=m*zg;
    Amat(2,2)=m-Zwdot;
    Amat(2,3)=-m*xg-Zqdot;
    Amat(3,1)=m*zg;
    Amat(3,2)=-m*xg-Mwdot;
    Amat(3,3)=Iyy-Mqdot;
    
    Bmat(1,1)=-Xuu*abs(up); % Matrix B
    Bmat(1,2)=m*qp-Xwq*qp;
    Bmat(1,3)=-m*xg*qp-Xqq*qp;
    Bmat(2,1)=-m*qp-Zuq*qp-Zuud*up*del;
    Bmat(2,2)=-Zww*abs(wp)-Zuw*up;
    Bmat(2,3)=-m*zg*qp-Zqq*abs(qp);
    Bmat(3,1)=-Muw*wp-Muud*up*del;
    Bmat(3,2)=-Mww*abs(wp);
    Bmat(3,3)=m*zg*wp+m*xg*up-Mqq*abs(qp)-Muq*up;
    
    exLoad(1,1)=Xprop;
    
    Cmat=Amat/dt+Bmat/2;    % Crank Nicholson Technique
    FF=(Amat/dt-Bmat/2)*solo+exLoad;
    
    soln=Cmat\FF;           % new solution
    
    diff=(soln(1)-solo(1))^2+(soln(1)-solo(1))^2+(soln(1)-solo(1))^2;
    ratio=sqrt(diff/(soln(1)^2+soln(2)^2+soln(3)^2));
    
    up=(solo(1)+soln(1))/2; 
    wp=(solo(2)+soln(2))/2; 
    qp=(solo(3)+soln(3))/2;    
    NoIt(itime+1)=it;
    ConvError(itime+1)=ratio;
    
    end  % iteration loop for convergence
    
    u(itime+1)=soln(1);
    w(itime+1)=soln(2);
    q(itime+1)=soln(3);
    
    x_position(itime+1)=x_position(itime)+(u(itime+1)+u(itime))*dt/2;
    z_position(itime+1)=z_position(itime)+(w(itime+1)+w(itime))*dt/2;
    q_position(itime+1)=q_position(itime)+(q(itime+1)+q(itime))*dt/2;
    
end  % end loop for itime


figure()
subplot(3,1,1)
plot(time,u)
xlabel('time (sec)')
ylabel('u (m/s)')
grid on; box on;
title('Velocity Plot')
subplot(3,1,2)
plot(time,w)
xlabel('time (sec)')
ylabel('w (m/s)')
grid on; box on;
subplot(3,1,3)
plot(time,q)
xlabel('time (sec)')
ylabel('q (rad/s)')
grid on; box on;

figure()
subplot(3,1,1)
plot(time,x_position)
xlabel('time (sec)')
ylabel('x position (m)')
grid on; box on;
title('Poistion Plot')
subplot(3,1,2)
plot(time,z_position)
xlabel('time (sec)')
ylabel('z position(m)')
grid on; box on;
subplot(3,1,3)
plot(time,q_position)
xlabel('time (sec)')
ylabel('pitch angle (ra)')
grid on; box on;

figure()
subplot(2,1,1)
plot(time,NoIt)
xlabel('time (sec)')
ylabel('No. of Iterations')
title('No. of Itertaions for Convergency')

subplot(2,1,2)
plot(time,ConvError)
xlabel('time (sec)')
ylabel('Ratio of Solutions')
title('Error at Convergence')
































