%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   2D Dynamic analysis of REMUS in xz-plane
%
%         19July2020
%
%  This program comuputes the dynamic motion of the REMUS underwater vehicles
%  subjected to various hydrodynamic loads in the xz-plane.
%    
%   [Amat]{acc}+[Bmat]{ve}={HydroLoad} is soleved using the
%   Explicit Euler technique
%   Here, {vel}={u; w; q} column vector
%         {acc}=(udot; wdot; qdot} column vector
%         {HydroLoad} is the hydrodynamic forces and moment
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
%   Fz: linear force along the z-axis
%   My: moment about the y-axis
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear('all');
close('all');

dt=0.0005;      % time step size (sec)
FinalTime=60;      % final time of simulation
del=20.0*(pi/180);       % fin angle
PlotVehicleLength = 0.75;

m=30;       % mass in kg
Iyy=3.45;   % mass moment of interia about y-axis in kg-m^2

xg=0;       % center of gravisty wrt center of bouyance in x-axis 
zg=1.96e-2; % center of gravisty wrt center of bouyance in z-axis (m)

Xuu=-1.62;          % Axial flow drag (kg/m)
Xudot=-9.30e-001;	% Added Mass (kg)
Xwq=-3.55e+001;     % Added Mass Cross-term (kg/rad)
Xqq=1.93e+000;      % Added Mass Cross-term (kg-m/rad)
Xprop=+3.86e+000;   % Propeller Thrust (N)
Zww=-1.31e+002;     % Cross-flow Drag (kg/m)
Zqq=-6.32e-001;     % Cross-flow Drag (kg-m/rad^2)
Zuw=-2.86e+001;  	% Body Lift Force and Fin Lift (kg/m)
Zwdot=-3.55e+001;   % Added Mass (kg)
Zqdot=-1.93e+000;   % Added Mass (kg-m/rad)
Zuq=-5.22e+000;     % Added Mass Cross-term and Fin Lift (kg/m)
Zuud=-9.64e+000;    % Fin Lift Force (kg/(m-rad))
Mww=+3.18e+000;     % Crossflow Drag (kg)
Mqq=-9.40e+000;     % Crossflow Drag (kg-m^2 / rad^2)
Muw=+2.40e+001;     % Body and Fin Lift and Munk Moment (kg)
Mwdot=-1.93e+000;   % Added mass (kg-m)
Mqdot=-4.88e+000;   % Added mass (kg-m^2 /rad)
Muq=-2.00e+000;     % Added Mass Cross Term and Fin Lift (kg-m/rad)
Muud=-6.15e+000;    % Fin Lift Moment (kg/rad)

ntime=FinalTime / dt;      % no of time steps

Amat=zeros(3,3);    % initialization of matrice Amat
Bmat=zeros(3,3);    % initialization of matrice Bmat
Hvec=zeros(3,1);  % initialization of column vector Hvec

u=zeros(ntime+1,1);     % velocity in x-axis
w=zeros(ntime+1,1);     % velocity in z-axis
q=zeros(ntime+1,1);     % angular velocity 

x_position=zeros(ntime+1,1);        % position in x-axis
z_position=zeros(ntime+1,1);        % position in z-axis
q_position=zeros(ntime+1,1);        % angular position

time=zeros(ntime+1,1);     % time

u(1)=0;     % inital x-velocity
w(1)=0;     % inital z-velocity
q(1)=0;     % inital angular velocity about y-axis

% Matrix A
Amat(1,1)=m-Xudot;
Amat(1,3)=m*zg;
Amat(2,2)=m-Zwdot;
Amat(2,3)=-m*xg-Zqdot;
Amat(3,1)=m*zg;
Amat(3,2)=-m*xg-Mwdot;
Amat(3,3)=Iyy-Mqdot;

for itime=1:ntime   % time increment loop
    
    itime
    time(itime+1)=time(itime)+dt;
    
    % Matrix B
    Bmat(1,1)=-Xuu*abs(u(itime));
    Bmat(1,2)=m*q(itime)-Xwq*q(itime);
    Bmat(1,3)=-m*xg*q(itime)-Xqq*q(itime);
    Bmat(2,1)=-m*q(itime)-Zuq*q(itime)-Zuud*u(itime)*del;
    Bmat(2,2)=-Zww*abs(w(itime))-Zuw*u(itime);
    Bmat(2,3)=-m*zg*q(itime)-Zqq*abs(q(itime));
    Bmat(3,1)=-Muw*w(itime)-Muud*u(itime)*del;
    Bmat(3,2)=-Mww*abs(w(itime));
    Bmat(3,3)=m*zg*w(itime)+m*xg*u(itime)-Mqq*abs(q(itime))-Muq*u(itime);
    
    Hvec(1,1)=Xprop;
    
    Solution = [u(itime); w(itime); q(itime)] + dt*( Amat \ (Hvec - Bmat*[u(itime); w(itime); q(itime)]) );
    u(itime+1) = Solution(1);
    w(itime+1) = Solution(2);
    q(itime+1) = Solution(3);
    
end  % end loop for itime

x_position = cumtrapz(time,u);
z_position = cumtrapz(time,w);
q_position = cumtrapz(time,q);

xBow = x_position + PlotVehicleLength/2*cos(q_position);
xStern = x_position - PlotVehicleLength/2*cos(q_position);
zBow = z_position + PlotVehicleLength/2*sin(q_position);
zStern = z_position - PlotVehicleLength/2*sin(q_position);

Length = sqrt( (xBow-xStern).^2 + (zBow-zStern).^2 );

figure()
subplot(3,1,1)
plot(time,u)
xlabel('time (sec)')
ylabel('u (m/s)')
grid on; box on;
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
print('Velocity_del_20deg','-dpdf','-r600','-bestfit');

figure()
subplot(3,1,1)
plot(time,x_position)
xlabel('time (sec)')
ylabel('x (m)')
grid on; box on;
subplot(3,1,2)
plot(time,z_position)
xlabel('time (sec)')
ylabel('z (m)')
grid on; box on;
subplot(3,1,3)
plot(time,(180/pi)*q_position)
xlabel('time (sec)')
ylabel('\theta (deg)')
grid on; box on;
print('Position_del_20deg','-dpdf','-r600','-bestfit');

figure()
hold on;
plot(x_position,z_position,'-')
for i = 1:1000:length(x_position)
    plot([xBow(i) xStern(i)],[zBow(i) zStern(i)],'-k')
end
xlabel('x (m)')
ylabel('z (m)')
axis equal
grid on; box on;
print('Trajectory_del_20deg','-dpdf','-r600','-bestfit');
