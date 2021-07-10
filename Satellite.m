function derivative = Satellite(t,state)

global BB invI I Bfieldmeasured pqrdotmeasured
global Bfieldnav pqrdotnav current

%%Earth and sat params
Earth
m = 2;


x = state(1);
y = state(2);
z = state(3);

q0123 = state(7:10);
pdot = state(11);
qdot = state(12);
rdot = state(13);
pqrdot = state(11:13);

%%translation
vel = state(4:6);
r = state(1:3);
rho = norm(r);
rhat = r/rho;
Fgrav = -(G*M*m/rho^2)*rhat;
acc = Fgrav/m;

%%Mag stuff
phiE = 0;
thetaE = acos(z/rho);
psiE = atan2(y,x);
latitude = 90-thetaE*180/pi;
longitude = psiE*180/pi;
rhokm = (rho)/1000;
[BN,BE,BD] = igrf('01-Jan-2020',latitude,longitude,rhokm,'geocentric');
%%%Convert NED (North East Down to X,Y,Z in ECI frame)
%%%First we need to create a rotation matrix from the NED frame to the 
%%%inertial frame
BNED = [BN;BE;BD]; 
BI = TIB(phiE,thetaE+pi,psiE)*BNED;
%BI = eye(3)*BNED;    
BB = TIBquat(q0123)'*BI;
%%%Convert to Tesla
BB = BB*1e-9;

%%Bfield and ang velocty measured from snsr
[Bfieldmeasured,pqrdotmeasured] = Sensor(BB,pqrdot);

%%Navigation Block (correction)
[Bfieldnav,pqrdotnav] = Navigation(Bfieldmeasured,pqrdotmeasured);


%%CONTROL BLOCK
current = Control(Bfieldnav,pqrdotnav);
Magtorquer_params
muB = current*n*A;      %%Mag field produced by torquer
LMN_magtorquers = cross(muB,BB);

%%%Rotational Dynamics
Inertia
PQRMAT = [0 -pdot -qdot -rdot;pdot 0 rdot -qdot;qdot -rdot 0 pdot;rdot qdot -pdot 0];
q0123dot = 0.5*PQRMAT*q0123;
H = I*pqrdot;
pqrddot = invI*(LMN_magtorquers - cross(pqrdot,H));

derivative = [vel;acc;q0123dot;pqrddot];

