
%%The following code and its dependencies have been created with the help
%%from Asstnt Professor of University of South Alabama, Carlos Montalvo.


%%Main
clear
clc
close all

disp('Sim started');

tic
global BB Bfieldmeasured pqrdotmeasured current
global Bfieldnav pqrdotnav Bfieldnavprev pqrdotnavprev

%%Earth and orbit params
Earth
Orbit
addpath './igrf'

%%Initial Params
x0 = Rad_mod;
y0 = 0;
z0 = 0;
xdot0 = 0;
ydot0 = Velocity*cos(Inclination);
zdot0 = Velocity*sin(Inclination);
p = 0;
q = 0;
r = 0;
pqr = [p q r];
q0123_0 = EulerAngles2Quaternions(pqr);
pdot = 0.175;
qdot = 0.175;
rdot = 0.175;
state = [x0 y0 z0 xdot0 ydot0 zdot0 q0123_0 pdot qdot rdot];


%%Time params
period = 2*pi*sqrt((Rad_mod^3)/mu);
no_of_orbs = 18;
tfinal = no_of_orbs*period;
timestep = 1;
tout = 0:timestep:tfinal;
stateout = zeros(length(tout),length(state));

%%Initialising mag matrices
Bxs = zeros(length(tout),1);    %%B field from IGRF
Bys = zeros(length(tout),1);
Bzs = zeros(length(tout),1);
Bxm = zeros(length(tout),1);    %%Polluted B field
Bym = zeros(length(tout),1);
Bzm = zeros(length(tout),1);
Bxn = zeros(length(tout),1);    %%Navigated B field
Byn = zeros(length(tout),1);
Bzn = zeros(length(tout),1);
Bfieldnavprev = [0;0;0];

%%Initialising polluted ang velocity matrices
pqrdot_measured = zeros(length(tout),3);
pqrdot_navigation = zeros(length(tout),3);
pqrdotnavprev = [0;0;0];

%%Initialising current matrix
currentmatrix = zeros(length(tout),3);


%%RK4
state = state';
for i = 1:length(tout)
    stateout(i,:) = state;
    k1 = Satellite(tout,state);
    k2 = Satellite(tout + timestep/2,state + k1*timestep/2);
    k3 = Satellite(tout + timestep/2,state + k2*timestep/2);
    k4 = Satellite(tout + timestep,state + k3*timestep);
    k = (1/6)*(k1 + 2*k2 + 2*k3 + k4);
    state = state + k*timestep;
    
    Bxs(i,1) = BB(1);
    Bys(i,1) = BB(2);
    Bzs(i,1) = BB(3);
    
    Bxm(i,1) = Bfieldmeasured(1);
    Bym(i,1) = Bfieldmeasured(2);
    Bzm(i,1) = Bfieldmeasured(3);
    
    Bxn(i,1) = Bfieldnav(1);
    Byn(i,1) = Bfieldnav(2);
    Bzn(i,1) = Bfieldnav(3);
    
    pqrdot_measured(i,1) = pqrdotmeasured(1);
    pqrdot_measured(i,2) = pqrdotmeasured(2);
    pqrdot_measured(i,3) = pqrdotmeasured(3);
    
    pqrdot_navigation(i,1) = pqrdotnav(1);
    pqrdot_navigation(i,2) = pqrdotnav(2);
    pqrdot_navigation(i,3) = pqrdotnav(3);
    
    currentmatrix(i,1) = current(1);
    currentmatrix(i,2) = current(2);
    currentmatrix(i,3) = current(3);
    
    if mod(i,10) == 0
        perc = (i/length(tout))*100;
        disp(['Percentage completed: ' num2str(perc)]);
    end
    
    
end

%%In kms
stateout(:,1:6) = stateout(:,1:6)/1000;

%%Extraction
pqrdot_out = stateout(:,11:13);
Bnorm = sqrt(Bxs.^2 + Bys.^2 +Bzs.^2);
Xout = stateout(:,1);
Yout = stateout(:,2);
Zout = stateout(:,3);

%%Display
disp('Simulation completed');
toc

%%plot x y z
fig0 = figure();
plot(tout,Xout,'r-');
hold on
plot(tout,Yout,'g-');
hold on
plot(tout,Zout,'b-');

%%plot Bx By Bz wrt BInertial frame
fig1 = figure();
plot(tout,Bxs,'b-','Linewidth',2);
hold on
plot(tout,Bys,'g-','LineWidth',2);
plot(tout,Bzs,'r-','LineWidth',2);
plot(tout,Bxm,'b--');
plot(tout,Bym,'g--');
plot(tout,Bzm,'r--');
plot(tout,Bxn,'k--');
plot(tout,Byn,'c--');
plot(tout,Bzn,'m--');
ylabel('B components in T');
xlabel('Time');

%%Plot Bnorm wrt Inertial frame
fig2 = figure();
plot(tout,Bnorm,'LineWidth',3);
xlabel('Time');
ylabel('B');
ylabel('Magnetic field in T');
xlabel('Time');

%%Plot Angular Velocity
fig3 = figure();
plot(tout,pqrdot_out,'LineWidth',2);
hold on 
%%plot(tout,pqrdot_measured,'--');
%%plot(tout,pqrdot_navigation,'--','LineWidth',2);
xlabel('Time');
ylabel('Angular Velocity');

%%Plot current
fig4 = figure();
plot(tout,currentmatrix);
