%% General Conditions %%

T0 = 200;               % Temperature (K)
Vn = 512;               % Atom size (m^3)
N = 30;                 % Number of Atoms
Vol = N * Vn;           % Total Volume (m^3)
side = Vol^(1.0/3.0);   % Length of Side of Simulation Volume (Angstrom)
dt = 0.01;             % Time Step (s)
MW = 0;                 % Molecular Weight (Kilograms/Mole)
Na = 6.022e+23;         % Avogadro's Number (Atoms/Mole)
eps = 136;              % Depth of the Potential Well (eV)              ***
sigma = 2.89e-10;       % Collision Diameter (m)
rCut = 15;              % Cut-off Distance (m)
maxStep = 3;            % Upper bound for iterations
U = 0;                  % Potential Energy (J)                          

%% Property Initialization %%
% Here we have matrices of N atoms/particles containing their mechanical
% propeties.

particles = Particle([]);

%% Initialize Positions %%

particles = PositionInitialization(N, side, particles);   % Get initial position for particles
particles = InitializeVelocity(particles, N, T0);
pos = zeros(N, 3);                                        % Pre-allocation of memory

for i = 1:3
   pos(:,i) = GetVectorProps(particles, N, i);            % Position vector for particles                
end

% Generate plot for initial position
plot3(pos(:,1), pos(:,2), pos(:,3), 'o', 'MarkerFaceColor', 'k');
xlabel('x'), ylabel('y'), zlabel('z'), title('Coordinates of Particles');
xlim([0 ceil(side)]), ylim([0 ceil(side)]), zlim([0 ceil(side)]);
grid on

%% Initial Computations and Correction Balancing Variables %%

for i = 1:N
   MW = MW + particles(i).SpecieData.elementWeight; 
end

%mass = (N / Na) * MW;                                      % (Kg)
mass = MW * 1e28 / Na / 1000;
rNbr = rCut + 3.0;
rNbr2 = rNbr * rNbr;


%% Neighbour List %%
% Compute distance between all pairs of particles

particles = ComputeDistance(particles, N, rNbr2);

%% Initialize Force %%

particles = ComputeForce(particles, N, 1);

%% Subroutines %%
% The following methods will compute thew new position of particles,
% evaluate the forces in the system, correct position computations based on
% repulsiong/attractive forces and gear corrector, apply the periodic
% boundary conditions, scale velocity, update neighbours list and print the
% results for each time step up until the max time.

fprintf('Time (s) Kinetic Energy (aJ) Lennard Jones Potential (aJ) Total Energy (aJ) Temperature (K)\n');
LJ = 0;
x = 0;
j = 1;

R = 8.314;                                                  % Gas Constant
Ek_Adjust = 1;                                              % Kinetic Energy Adjustment Factor
velocity_scale = 1;                                         % Constant used to maintain constant T
bounds = zeros(1, 2); bounds(1) = 100; bounds(2) = 0.5;

for t = 1:dt:maxStep
   particles = PositionPredictor(particles, N, dt, velocity_scale);
   particles = BoundaryCondition(N, particles, side);
   particles = ComputeForce(particles, N, velocity_scale);
   particles = LennardJonesEvaluator(particles, N, sigma, eps);
   particles = ComputeDistance(particles, N, rNbr2);
   
   ULJ = TotalLennardJones(N, particles);
   TotVelocity = TotalVelocity(N, particles);
   
   KE = (1 / 2) * mass * TotVelocity^2;
   T = ComputeSystemTemperature(KE, Ek_Adjust, Na, R);
   velocity_scale = TemperatureScale(T, T0, bounds);
   TotE = KE + ULJ;
   fprintf('%2.3f %15.3f %25.3f %23.3f %10.3f\n', t - 1, KE, ULJ, TotE, T);
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