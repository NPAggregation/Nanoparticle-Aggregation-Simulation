%% General Conditions %%
T = 300;                % Temperature (K)
Vn = 512;               % Atom size (Angstroms Cubed/Atom)
N = 30;                 % Number of Atoms
Vol = N * Vn;           % Total Volume (Angstroms^3)
side = Vol^(1.0/3.0);   % Length of Side of Simulation Volume (Angstrom)
dt = 0.3;               % Time Step (fs)
MW = 1.0;               % Molecular Weight (Grams/mole)
halfSide = 0.5 * side;  % Half of the Side (Angstrom)
density = 1 / Vn;       % Molar Density (Atom/Angstrom)
Na = 6.022e+23;         % Avogadro's Number (Atoms)
eps = 136;              % Depth of the Potential Well (eV)              ***
sigma = 2.89;           % Collision Diameter (Angstroms)                
rCut = 15;              % Cut-off Distance (Angstroms)
Q = 0;                  % Charge of Atom (Coulomb)                      ***
e0 = 8.854187;          % Vacuum Permittivity (Farads/metre)            ***
kb = 1.38066e-5;        % Boltzmann's Constant (aJ/molecule/K)          ***
maxStep = 50;           % Upper bound for iterations
%U = 0;                 % Potential Energy (J)                          ***
sampleIntval = 5;       % Sampling Interval
nbrIntval = 10;         % Neighbor's List Update Interval
writeIntval = 100;      % Writing Interval
nProp = 4;              % Number of Properties

%% Property Initialization %%
% Here we have matrices of N atoms/particles containing their mechanical
% propeties.

particles = Particle([]); 

%% Initial Computations and Correction Balancing Variables %%
mass = MW / Na / 1000 * 1e28;               % (kg / molecule)
Tcorr = (3.0 * N * kb * T) / mass;          % Temperature Correction Factor (Angstrom/fs)^2
dtva = [dt dt^2 dt^3 dt^4 dt^5];            % Gear Correction Factor
fv = [1 2 6 24 120];                        % Vector of Factorials
dtv = dtva./fv;                             % Gear Correction Factor
rNbr = rCut + 3.0;                          % Radius Defining if Particle is a Neighbours (Angstrom)
rCut2 = rCut * rCut;
rNbr2 = rNbr * rNbr;
rCut3 = rCut^3;
rCut9 = rCut^9;

% Corrector Coefficients for Gear Using Dimensioned Variables
gear = [(3.0 / 2.0) (251.0 / 360.0) (1.0) (11.0 / 18.0) (1.0 / 6.0) (1.0 / 60.0)];
dtv6 = [1 dt dt^2 dt^3 dt^4 dt^5];
fv6 = [1 1 2 6 24 120];
alpha(1:6) = gear(1:6)./(dtv6(1:6).*fv6);
alpha = (alpha * dt^2) / 2;

%% Initialize Positions %%
particles = PositionInitialization(N, side, particles);   % Get initial position for particles
pos = zeros(N, 3);                                        % Pre-allocation of memory
for i = 1:3
   pos(:,i) = GetVectorProps(particles, N, i);            % Position vector for particles                
end

% Generate plot for initial position
plot3(pos(:,1), pos(:,2), pos(:,3), 'o', 'MarkerFaceColor', 'k');
xlabel('x'), ylabel('y'), zlabel('z'), title('Coordinates of Particles');
xlim([0 ceil(side)]), ylim([0 ceil(side)]), zlim([0 ceil(side)]);
grid on

%% Initial Velocity %%
% Initialize the velocity of particles with zerno net momentum in the
% system.

particles = VelocityInitialize(particles, N, Tcorr);

%% Neighbour List %%
% Compute distance between all pairs of particles

particles = ComputeDistance(particles, N, rNbr2);


%% Subroutines %%
% The following methods will compute thew new position of particles,
% evaluate the forces in the system, correct position computations based on
% repulsiong/attractive forces and gear corrector, apply the periodic
% boundary conditions, scale velocity, update neighbours list and print the
% results for each time step up until the max time.
fprintf('Time (s) Kinetic Energy (aJ) Lennard Jones Potential (aJ) Total Energy (aJ)\n');
U = 0;
LJ = 0;
x = 0;
j = 1;
for t = 1:dt:maxStep
   particles = PositionPredictor(particles, N, dtv);
   particles = VelocityScale(N, particles, Tcorr);
   particles = BoundaryCondition(N, particles, side);
   particles = ForceEvaluator(particles, N, sigma, eps);
   particles = LennardJonesEvaluator(particles, N, sigma, eps);
   if (mod(t, 1) == 1)
       particles = RepelParticles(particles, N);
   end
   ULJ = TotalLennardJones(N, particles);
   particles = ComputeDistance(particles, N, rNbr2);
   
   TotVelocity = TotalVelocity(N, particles);
   KE = (1 / 2) * mass * TotVelocity^2;
   TotE = KE + ULJ; 
   fprintf('%2.1f %15.3f %25.3f %23.3f\n', t, KE, ULJ, TotE);
   U(j) = ULJ;
   LJ(j) = particles(1).LJ(2);
   x(j) = particles(1).NeighborList(2);
   j = j + 1;
   
   pos = zeros(N, 3);
   for i = 1:3
       pos(:, i) = GetVectorProps(particles, N, i);
   end
   
   pause(0.1)
   plot3(pos(:,1), pos(:,2), pos(:,3), 'o', 'MarkerFaceColor', 'k');
   xlabel('x'), ylabel('y'), zlabel('z'), title('Coordinates of Particles');
   xlim([0 ceil(side)]), ylim([0 ceil(side)]), zlim([0 ceil(side)]);
   grid on
end
t = 1:dt:maxStep;
figure
plot(t, U); title('Total Lennard Jones vs Time')
figure
plot(t, LJ); title('Lennard Jones Between i, j Pair vs Time')
figure
plot(x, LJ); title('Lennard Jones vs Distance')