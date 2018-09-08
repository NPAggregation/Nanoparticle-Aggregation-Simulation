clear, clc
%% Last edit: 29/01/2018
%% GENERAL CONDITIONS %%
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
maxstp = 2000;
U = 0;
ksamp = 1;      % sampling interval
knbr = 10;      % neighbor list update interval
kwrite = 100;   % writing interval
nprop = 4; %Properties
%% PROPERTY INITIALIZATION %%
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

%% COMPUTATION AND CORRECTIONS %%
mass = MW/Nav/1000*1.0e+28; % (1e-28*kg/molecule)
tfac = 3.0*N*kb*T/mass; % Temperature correction factor (Angstrom/fs)^2
dtva = [dt dt*dt dt*dt*dt dt*dt*dt*dt dt*dt*dt*dt*dt];
fv = [1 2 6 24 120]; % vector of factorials
dtv = dtva./fv;
rnbr = rcut + 3.0;
rcut2 = rcut*rcut;
rnbr2 = rnbr*rnbr;
rcut3 = rcut^3;
rcut9 = rcut^9;

% corrector coefficients for Gear using dimensioned variables
gear = [3.0/20.0 251.0/360.0 1.0 11.0/18.0 1.0/6.0 1.0/60.0];
dtv6 = [1 dt dt*dt dt*dt*dt dt*dt*dt*dt dt*dt*dt*dt*dt];
fv6 = [1 1 2 6 24 120]; % vector of factorials
alpha(1:6) = gear(1:6)./dtv6(1:6).*fv6;
alpha = alpha*dt^2/2;

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
plot3(r(:,1),r(:,2),r(:,3),'or')
ylabel('y'), xlabel('x'),zlabel('z');
ylim([0,ceil(side)]),xlim([0,ceil(side)]),zlim([0,ceil(side)]);
set(gca,'Color',[1 1 1])
grid on
title('Display of r')
%% Initial Potential Calculator

ULJ = 4*eps*((sigma/rcut)^12-(sigma/rcut)^6); % Lennard-Jones potential (eV)
Ucoulomb = Q^2/(4*pi*e0*rcut); % Coulomb potential
%Bonding Potential
Utot = ULJ+Ucoulomb; %Total Energy (other terms to be added)


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

%% Neighbor List
Nnbr = 0;
for i = 1:1:N
    for j = i+1:1:N
        dis(1:3) = r(i,1:3) - r(j,1:3); % STORES DISTANCE BETWEEN PAIR (i,j)
        for k = 1:1:3
            if (dis(k) > sideh); dis(k) = dis(k) - side; end;
            if (dis(k) < -sideh); dis(k) = dis(k) + side; end;
        end
        dis2 = sum(dis.*dis);
        if (dis2 <= rnbr2)
            Nnbr = Nnbr + 1;
            Nnbrlist(Nnbr,1) = i;
            Nnbrlist(Nnbr,2) = j;
        end
    end
end
if (Nnbr == 0)
    Nnbrlist = zeros(1,1);
end
%% Subroutines
for istep = 1:1:maxstp
    %% New Position Predictor
    r = r + v *dtv(1) + dtv(2)*a + dtv(3)*d3 + dtv(4)*d4 + dtv(5)*d5;
    v = v + a *dtv(1) + dtv(2)*d3 + dtv(3)*d4 + dtv(4)*d5;
    a = a + d3*dtv(1) + dtv(2)*d4 + dtv(3)*d5;
    d3 = d3 + d4*dtv(1) + dtv(2)*d5;
    d4 = d4 + d5*dtv(1);
    
    %% Evaluate Forces
    sig6 = sigma^6; sig12 = sigma^12;
    U=0; f = zeros(N,3); % reinitialize every time
    for n = 1:1:Nnbr
        i = Nnbrlist(n,1);
        j = Nnbrlist(n,2);
        dis(1:3) = r(i,1:3) - r(j,1:3); % STORES DISTANCE BETWEEN PAIR (i,j)
        for k = 1:1:3
            if (dis(k) > sideh); dis(k) = dis(k) - side; end;
            if (dis(k) < -sideh); dis(k) = dis(k) + side; end;
        end
        dis2 = sum(dis.*dis);
        if (dis2 <= rcut2) % IF THE CUT OFF CRITERIA IS SATISFIED
            dis2i = 1.0/dis2;
            dis6i = dis2i*dis2i*dis2i;
            dis12i = dis6i*dis6i;
            U = U + ( sig12*dis12i - sig6*dis6i );
            fterm = (2.0*sig12*dis12i - sig6*dis6i )*dis2i;
            %fprintf(1, 'n %i i %i j %i dis2 %e fterm %e %e %e\n',n,i,j,dis2,fterm);
            f(i,1:3) = f(i,1:3) + fterm.*dis(1:3);
            f(j,1:3) = f(j,1:3) - fterm.*dis(1:3);
        end
    end
    f = f*24.0*eps;
    U = U*4.0*eps;
    if istep == 1
        a = f/mass;
    end
    
    %% Position Corrector (IN COMMENTS BECAUSE NOT CURRENTLY WORKING)
%         for i = 1:1:N
%             errvec(1:3) = ( f(i,1:3)/mass - a(i,1:3) );
%             r(i,1:3) = r(i,1:3) + errvec(1:3)*alpha(1);
%             v(i,1:3) = v(i,1:3) + errvec(1:3)*alpha(2);
%             a(i,1:3) = a(i,1:3) + errvec(1:3)*alpha(3);
%             d3(i,1:3) = d3(i,1:3) + errvec(1:3)*alpha(4);
%             d4(i,1:3) = d4(i,1:3) + errvec(1:3)*alpha(5);
%             d5(i,1:3) = d5(i,1:3) + errvec(1:3)*alpha(6);
%         end
    %% Periodic Boundary Conditions
    for i = 1:1:N
        for j = 1:1:3
            if (r(i,j) > side); r(i,j) = r(i,j) - side; end;
            if (r(i,j) < 0.0); r(i,j) = r(i,j) + side; end;
        end
    end
    %% Velocity Scaling
    sumvsq = sum(sum(v.*v,1) );
    fac = sqrt(tfac/sumvsq);
    v = v*fac;
    %% Neighbor List
    Nnbr = 0;
    for i = 1:1:N
        for j = i+1:1:N
            dis(1:3) = r(i,1:3) - r(j,1:3);
            for k = 1:1:3
                if (dis(k) > sideh); dis(k) = dis(k) - side; end;
                if (dis(k) < -sideh); dis(k) = dis(k) + side; end;
            end
            dis2 = sum(dis.*dis);
            if (dis2 <= rnbr2)
                Nnbr = Nnbr + 1;
                Nnbrlist(Nnbr,1) = i;
                Nnbrlist(Nnbr,2) = j;
            end
        end
    end
    if (Nnbr == 0)
        Nnbrlist = zeros(1,1);
    end
    %% Property Printer
    props = zeros(nprop,6);
    den = floor(maxstp/ksamp);
    props(1:nprop,4) = props(1:nprop,2)/den;
    propname{1} = 'Kinetic Energy (aJ) ';
    propname{2} = 'Potential Energy (aJ) ';
    propname{3} = 'Total Energy (aJ) '; %%Plot the energy as a function of time
    propname{4} = 'Temperature (K) ';
    sumvsq = sum(sum(v.*v,1) );
    KE = 0.5*mass*sumvsq; % (aJ)
    T = 2/(3*N*kb)*KE;
    %
    % get instantaneous values
    %
    % property 1: total kinetic energy
    % property 2: total potential energy
    % property 3: total energy
    % property 4: temperature
    props(1,1) = KE;
    props(2,1) = U;
    props(3,1) = KE+U;
    props(4,1) = T;
    %% Plot & Displays
    % Plotting the r matrix
    if (mod(istep,kwrite) == 0)
        figure
        plot3(r(:,1),r(:,2),r(:,3),'or')
        ylabel('y'), xlabel('x'),zlabel('z');
        ylim([0,ceil(side)]),xlim([0,ceil(side)]),zlim([0,ceil(side)]);
        set(gca,'Color',[1 1 1])
        grid on
        title('Display of r')
        fprintf(1,'istep %i K %.5f \tU %.2f \tTOT %.2f \tT %.2f \n',istep,props(1:4,1));
    end
    
end