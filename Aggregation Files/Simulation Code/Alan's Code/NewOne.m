clear, clc
%% Last edit: 22/01/2018
%%%%%%%%%%%%%%%%GENERAL CONDITIONS%%%%%%%%%%%%%%%%%%%%%
T = 300;    % Temperature (K)
Vn = 512.0;     % Atom size, in Angstroms cubed / atom
N = 27;     % Number of atoms
Vol = N*Vn;     % total volume (Angstroms^3)
side = Vol^(1.0/3.0);   % length of side of simulation volume (Angstrom)
dt = 1; %Time Step (fs)
MW = 16.0;      % molecular weight (grams/mole)
sideh = 0.5*side;   % half of the side
density = 1/Vn;     % molar density
Nav = 6.022e+23; % Avogadro's Number
eps = 136; % Depth of the potential well (K or eV) to be decided
sigma = 3.884; % Collision Diameter (Angstroms)
rcut = 15; % Cut-off Distance
Q = 0; % Charge of the atom (Coulomb)
e0 = 8.854187; % Vacuum permittivity (farads per metre)
kb = 1.38066e-5; % Boltzmann's constant (aJ/molecule/K)

%%%%%%%%%%%%%%%PROPERTY INITIALIZATION%%%%%%%%%%%%%%%%%
% (x,y,z) coordinates 
% N by 3 matrices
r = zeros(N,3); % position
v = zeros(N,3); % velocity
a = zeros(N,3); % acceleration
d3 = zeros(N,3); % third derivative
d4 = zeros(N,3); % fourth derivative
d5 = zeros(N,3); % fifth derivative
f = zeros(N,3); % force 
p = zeros(N,3); % momentum

%%%%%%%%%%%%%COMPUTATION AND CORRECTIONS%%%%%%%%%%%%%%%%
mass = MW/Nav/1000*1.0e+28; % (1e-28*kg/molecule)
tfac = 3.0*N*kb*T/mass; % Temperature correction factor (Angstrom/fs)^2 
dtva = [dt dt*dt dt*dt*dt dt*dt*dt*dt dt*dt*dt*dt*dt];
fv = [1 2 6 24 120]; % vector of factorials
dtv = dtva./fv;
%% Initial Positions

ni = ceil(N^(1.0/3.0));
ncount = 0;
dx = side/ni;
for ix = 1:1:ni
    for iy = 1:1:ni
        for iz = 1:1:ni
            if (ncount <= N)
                ncount = ncount + 1;
                r(ncount,1) = dx*ix;
                r(ncount,2) = dx*iy;
                r(ncount,3) = dx*iz;
            end
        end
    end
end

%% Initial Potential Calculator

ULJ = 4*eps*((sigma/rcut)^12-(sigma/rcut)^6); % Lennard-Jones potential (eV)
Ucoulomb = Q^2/(4*pi*e0*rcut); % Coulomb potential
%Bonding Potential
U = ULJ+Ucoulomb; %Total Energy (other terms to be added)


%% Initial Velocity


v = rand(N,3); %Assign Random velocities (Between 0 and 1)
v = 2*v - 1; %Make it between -1 and 1 
% enforce zero net momentum
for i = 1:1:3
    sumv(i) = sum(v(1:N,i));
    v(1:N,i) = v(1:N,i) - sumv(i)/N;
end
% scale initial velocities to set point temperature
sumvsq = sum(sum(v.*v,1) );
fac = sqrt(tfac/sumvsq);
v = v*fac;

%% New Position Predictor

r = r + v *dtv(1) + dtv(2)*a + dtv(3)*d3 + dtv(4)*d4 + dtv(5)*d5;
v = v + a *dtv(1) + dtv(2)*d3 + dtv(3)*d4 + dtv(4)*d5;
a = a + d3*dtv(1) + dtv(2)*d4 + dtv(3)*d5;
d3 = d3 + d4*dtv(1) + dtv(2)*d5;
d4 = d4 + d5*dtv(1);