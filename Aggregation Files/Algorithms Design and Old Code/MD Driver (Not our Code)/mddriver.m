function mddriver
%
% This code performs molecular dynamics simulations
% in the canonical ensemble (specify T, V, and N)
%
% Author: David Keffer
% Department of Chemical Engineering, University of TN
% Last Updated: September 19, 2001
%
%global maxstp kmsd N dt

% Last edit from Alan Daou: June 25, 2017
%************************************************************
% PROGRAM INITIALIZATION
%************************************************************
%
% This code uses length units of Angstroms (1.0e-10 s)
% time = fs (1.0e-15 s)
% mass = (1.0e-28 kg)
% energy = aJ (1.0e-18 J)
% Temperature = K
%
%
% Specify thermodynamic state
%
T = 300;    % Temperature (K)
Vn = 512.0;     % Angstroms cubed / molecule
N = 27;     % Number of molecules

%
% Specify Numerical Algorithm Parameters
%
maxeqb = 2500;      % Number of time steps during equilibration
maxstp = 2000;      % Number of time steps during data production
dt = 1.0e-0;        % size of time step (fs)

%
%   Specify pairwise potential parameters
%

sig = 3.884;    % collision diameter (Angstroms)
eps = 137;      % well depth (K)
MW = 16.0;      % molecular weight (grams/mole)
rcut = 15;      % cut-off distance for potential (Angstroms)
%
% Specify sampling intervals
%
nprop = 7;      % number of properties
ksamp = 1;      % sampling interval
knbr = 10;      % neighbor list update interval
kwrite = 100;   % writing interval
kmsd = 10;      % position save for mean square displacement
rnbr = rcut + 3.0;
fid_msd = fopen('md_msd.out','w');

%
% props
% first index is property
% property 1: total kinetic energy
% property 2: total potential energy
% property 3: total energy
% property 4: temperature
% property 5: total x-momentum
% property 6: total y-momentum
% property 7: total z-momentum
%
% second index is
% 1: instantaneous value
% 2: sum
% 3: sum of squares
% 4: average
% 5: variance
% 6: standard deviation
%
props = zeros(nprop,6);
%
% Initialize vectors
%
% first index of r is over molecules
% second index of r is over dimensionality (x,y,z)
r = zeros(N,3); % position
v = zeros(N,3); % velocity
a = zeros(N,3); % acceleration
d3 = zeros(N,3); % third derivative
d4 = zeros(N,3); % fourth derivative
d5 = zeros(N,3); % fifth derivative
f = zeros(N,3); % force
rwopbc = zeros(N,3); % position w/o pbc

%
% additional simulation parameters
%
lmsd = 1;   % logical variable for mean square displacement
lscale = 1;     % logical variable for temperature scaling
%************************************************************
% INITIALIZATION PART TWO
%************************************************************
%
% compute a few parameters
%
dt2 = dt*dt;
dt2h = 0.5*dt2;
Vol = N*Vn;     % total volume (Angstroms^3)
side = Vol^(1.0/3.0);   % length of side of simulation volume (Angstrom)
sideh = 0.5*side;   % half of the side
density = 1/Vn;     % molar density
sig6 = sig^6;
sig12 = sig^12;
rcut2 = rcut*rcut;
rnbr2 = rnbr*rnbr;
% stuff for long range energy correction
rcut3 = rcut^3;
rcut9 = rcut^9;
kb = 1.38066e-5; % Boltzmann's constant (aJ/molecule/K) CHECK
eps = eps*kb;
Ulong = N*8*eps*pi*density*( sig12/(9.0*rcut9) - sig6/(3.0*rcut3) ); % CHECK
% temperature factor for velocity scaling
Nav = 6.022e+23; % Avogadro's Number
mass = MW/Nav/1000*1.0e+28; % (1e-28*kg/molecule)
tfac = 3.0*N*kb*T/mass; % (Angstrom/fs)^2
% correction factors for numerical algorithm
fv = [1 2 6 24 120]; % vector of factorials
dtva = [dt dt*dt dt*dt*dt dt*dt*dt*dt dt*dt*dt*dt*dt];
dtv = dtva./fv;
% corrector coefficients for Gear using dimensioned variables
gear = [3.0/20.0 251.0/360.0 1.0 11.0/18.0 1.0/6.0 1.0/60.0];
dtv6 = [1 dt dt*dt dt*dt*dt dt*dt*dt*dt dt*dt*dt*dt*dt];
fv6 = [1 1 2 6 24 120]; % vector of factorials
alpha(1:6) = gear(1:6)./dtv6(1:6).*fv6;
alpha = alpha*dt^2/2;
%
% assign initial positions of molecules in FCC crystal structure
%
[r,rwopbc] = funk_ipos(N,side,r,rwopbc);
%
% assign initial velocities
%
v = funk_ivel(N,v,T,tfac);
%
% create neighbor list
%
[Nnbr,Nnbrlist] = funk_mknbr(N,r,rnbr2,side,sideh);
%
% evaluate initial forces and potential energy
%
[f,U] = funk_force(N,r,rcut2,side,sideh,Nnbr,Nnbrlist,sig6,sig12,eps);
a = f/mass; % initial acceleration
[props] = funk_getprops(N,v,mass,T,kb,U,props,nprop);
fprintf(1,'istep %i K %e U %e TOT %e T %e \n',0,props(1:4,1));
% display(r)
% figure
% Plotting the r matrix
plot3(r(:,1),r(:,2),r(:,3),'or')
ylabel('y'), xlabel('x'),zlabel('z');
ylim([0,ceil(side)]),xlim([0,ceil(side)]),zlim([0,ceil(side)]);
set(gca,'Color',[1 1 1])
grid on
title('Display of r')
F(1)=getframe(gcf);
props = zeros(nprop,6);
%************************************************************
% EQUILIBRATION
%************************************************************
for istep = 1:1:maxeqb
    % predict new positions
    [r,v,a,d3,d4,rwopbc] = predictor(N,r,rwopbc,v,a,d3,d4,d5,dtv);
    % evaluate forces and potential energy
    [f,U] = funk_force(N,r,rcut2,side,sideh,Nnbr,Nnbrlist,sig6,sig12,eps);
    % correct new positions
    [r,v,a,d3,d4,d5,rwopbc] = corrector(N,r,rwopbc,v,a,d3,d4,d5,f,dt2h,alpha,mass);
    %a = f/mass;
    % apply periodic boundary conditions
    r = pbc(N,r,side);
    % scale velocities
    if (lscale == 1)
        v = funk_scalev(N,v,T,tfac);
    end
    % update neighbor list
    if (mod(istep,knbr) == 0)
        [Nnbr,Nnbrlist] = funk_mknbr(N,r,rnbr2,side,sideh);
    end
    % sample properties
    if (mod(istep,ksamp) == 0)
        [props] = funk_getprops(N,v,mass,T,kb,U,props,nprop);
    end
    % save positions for mean square displacement
    if (mod(istep,kwrite) == 0)
        fprintf(1,'istep %i K %e U %e TOT %e T %e \n',istep,props(1:4,1));
        %         figure
        %         display(r)
        % Plotting the r matrix
        % Replace r by rwopbc if wanted
        plot3(r(:,1),r(:,2),r(:,3),'or')
        ylabel('y'), xlabel('x'),zlabel('z');
        ylim([0,ceil(side)]),xlim([0,ceil(side)]),zlim([0,ceil(side)]);
        set(gca,'Color',[1 1 1])
        grid on
        title('Display of r (equilibration)')
        for i=2:1:26
            F(i)=getframe(gcf); % save each figure produced
        end
    end
end
movie(F) %plays saved figures as movie
figure
%
% write equilibration results
%
if (maxeqb > ksamp)
    [props] = funk_report(N,props,nprop,maxeqb,ksamp);
end
%************************************************************
% PRODUCTION
%************************************************************
props = zeros(nprop,6);
lscale = 0;
if (lmsd)
    funk_msd(N,rwopbc,fid_msd);
end

for istep = 1:1:maxstp
    % predict new positions
    [r,v,a,d3,d4,rwopbc] = predictor(N,r,rwopbc,v,a,d3,d4,d5,dtv);
    % evaluate forces and potential energy
    [f,U] = funk_force(N,r,rcut2,side,sideh,Nnbr,Nnbrlist,sig6,sig12,eps);
    % correct new positions
    [r,v,a,d3,d4,d5,rwopbc] = corrector(N,r,rwopbc,v,a,d3,d4,d5,f,dt2h,alpha,mass);
    % % apply periodic boundary conditions
    r = pbc(N,r,side);
    % scale velocities
    if (lscale == 1)
        % v = scalev();
    end
    % update neighbor list
    if (mod(istep,knbr) == 0)
        [Nnbr,Nnbrlist] = funk_mknbr(N,r,rnbr2,side,sideh);
    end
    % sample properties
    if (mod(istep,ksamp) == 0)
        [props] = funk_getprops(N,v,mass,T,kb,U,props,nprop);
    end
    % save positions for mean square displacement
    if (mod(istep,kwrite) == 0)
        fprintf(1,'istep %i K %e U %e TOT %e T %e \n',istep,props(1:4,1));
        % display(r) if need is
        % Plotting the r matrix
        % Replace r by rwopbc if wanted
        plot3(r(:,1),r(:,2),r(:,3),'or')
        ylabel('y'), xlabel('x'),zlabel('z');
        ylim([0,ceil(side)]),xlim([0,ceil(side)]),zlim([0,ceil(side)]);
        set(gca,'Color',[1 1 1])
        grid on
        title('Display of r (production)')
        for i=1:1:20
            J(i)=getframe(gcf); % save each figure produced
        end
    end
    if (lmsd)
        if (mod(istep,kmsd) == 0)
            funk_msd(N,rwopbc,fid_msd);
        end
    end
end
movie(J) %plays saved figures as movie
%
if (maxstp > ksamp)
    [props] = funk_report(N,props,nprop,maxstp,ksamp);
end
fclose(fid_msd);



%************************************************************
% SUBROUTINES
%************************************************************

%
% funk_ipos: assigns initial positions
%
function [r,rwopbc] = funk_ipos(N,side,r,rwopbc);
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

rwopbc = r;
%
% funk_ivel: assigns initial velocities
%
function v = funk_ivel(N,v,T,tfac);
%
v = rand(N,3); % random velocities from 0 to 1
v = 2*v - 1.0; % random velocities from -1 to 1
% enforce zero net momentum
for i = 1:1:3
    sumv(i) = sum(v(1:N,i));
    v(1:N,i) = v(1:N,i) - sumv(i)/N;
end
% scale initial velocities to set point temperature
sumvsq = sum(sum(v.*v,1) );
fac = sqrt(tfac/sumvsq);
v = v*fac;
%
% funk_mknbr: create neighbor list
%
function [Nnbr,Nnbrlist] = funk_mknbr(N,r,rnbr2,side,sideh);
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

%
% funk_force: evaluate forces
%

function [f,U] = funk_force(N,r,rcut2,side,sideh,Nnbr,Nnbrlist,sig6,sig12,eps);
f = zeros(N,3); % force
U = 0.0; % potential energy
for n = 1:1:Nnbr
    i = Nnbrlist(n,1);
    j = Nnbrlist(n,2);
    dis(1:3) = r(i,1:3) - r(j,1:3);
    for k = 1:1:3
        if (dis(k) > sideh); dis(k) = dis(k) - side; end;
        if (dis(k) < -sideh); dis(k) = dis(k) + side; end;
    end
    dis2 = sum(dis.*dis);
    if (dis2 <= rcut2)
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

%
% predict new positions
%
function [r,v,a,d3,d4,rwopbc] = predictor(N,r,rwopbc,v,a,d3,d4,d5,dtv)
rwopbc = rwopbc + v *dtv(1) + dtv(2)*a + dtv(3)*d3 + dtv(4)*d4 + dtv(5)*d5;
r = r + v *dtv(1) + dtv(2)*a + dtv(3)*d3 + dtv(4)*d4 + dtv(5)*d5;
v = v + a *dtv(1) + dtv(2)*d3 + dtv(3)*d4 + dtv(4)*d5;
a = a + d3*dtv(1) + dtv(2)*d4 + dtv(3)*d5;
d3 = d3 + d4*dtv(1) + dtv(2)*d5;
d4 = d4 + d5*dtv(1);

%
% correct new positions
%
function [r,v,a,d3,d4,d5,rwopbc] = corrector(N,r,rwopbc,v,a,d3,d4,d5,f,dt2h,alpha,mass);
for i = 1:1:N
    errvec(1:3) = ( f(i,1:3)/mass - a(i,1:3) );
    rwopbc(i,1:3) = rwopbc(i,1:3) + errvec(1:3)*alpha(1);
    r(i,1:3) = r(i,1:3) + errvec(1:3)*alpha(1);
    v(i,1:3) = v(i,1:3) + errvec(1:3)*alpha(2);
    a(i,1:3) = a(i,1:3) + errvec(1:3)*alpha(3);
    d3(i,1:3) = d3(i,1:3) + errvec(1:3)*alpha(4);
    d4(i,1:3) = d4(i,1:3) + errvec(1:3)*alpha(5);
    d5(i,1:3) = d5(i,1:3) + errvec(1:3)*alpha(6);
end
%
% apply periodic boundary conditions
%
function r = pbc(N,r,side);
for i = 1:1:N
    for j = 1:1:3
        if (r(i,j) > side); r(i,j) = r(i,j) - side; end;
        if (r(i,j) < 0.0); r(i,j) = r(i,j) + side; end;
    end
end

%
% funk_scalev: scale velocities
%
function v = funk_scalev(N,v,T,tfac);
% scale velocities to set point temperature
sumvsq = sum(sum(v.*v,1) );
fac = sqrt(tfac/sumvsq);
v = v*fac;


%
% calculate properties for sampling

%
function [props] = funk_getprops(N,v,mass,T,kb,U, props,nprop);
% props
% first index is property
% property 1: total kinetic energy
% property 2: total potential energy
% property 3: total energy
% property 4: temperature
% property 5: total x-momentum
% property 6: total y-momentum
% property 7: total z-momentum
%
% second index is
% 1: instantaneous value
% 2: sum
% 3: sum of squares
% 4: average
% 5: variance
% 6: standard deviation
%
sumvsq = sum(sum(v.*v,1) );
KE = 0.5*mass*sumvsq; % (aJ)
T = 2/(3*N*kb)*KE;
%
% get instantaneous values
%
props(1,1) = KE;
props(2,1) = U;
props(3,1) = KE+U;
props(4,1) = T;
props(5,1) = mass*sum(v(:,1));
props(6,1) = mass*sum(v(:,2));
props(7,1) = mass*sum(v(:,3));
props(1:nprop,2) = props(1:nprop,2) + props(1:nprop,1);
props(1:nprop,3) = props(1:nprop,3) + props(1:nprop,1).*props(1:nprop,1);
%
% calculate and report simulation statistics
%
function [props] = funk_report(N,props,nprop,maxeqb,ksamp);
den = floor(maxeqb/ksamp);
props(1:nprop,4) = props(1:nprop,2)/den;
props(1:nprop,5) = props(1:nprop,3)/den - props(1:nprop,4).^2;
props(1:nprop,6) = sqrt(props(1:nprop,5));
propname{1} = 'Kinetic Energy (aJ) ';
propname{2} = 'Potential Energy (aJ) ';
propname{3} = 'Total Energy (aJ) '; %%Plot the energy as a function of time
propname{4} = 'Temperature (K) ';
propname{5} = 'x-Momentum ';
propname{6} = 'y-Momentum ';
propname{7} = 'z-Momentum ';
fprintf(1,' property instant average standard deviation \n');
for i = 1:1:nprop
    fprintf(1,' %s %e %e %e \n', propname{i}, props(i,1),props(i,4),props(i,6))
end

%
% save positions for mean square displacement calculations
%
function funk_msd(N,rwopbc,fid_msd);
for i = 1:1:N
    fprintf(fid_msd,' %e %e %e \n',rwopbc(i,1:3));
end