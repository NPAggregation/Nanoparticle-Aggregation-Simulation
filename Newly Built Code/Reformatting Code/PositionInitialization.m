%% Function to Initialize Initial Positions of Particles %%
% Function intializes position of each particle

function posMatrix = PositionInitialization(N, side, pos)
n = ceil(N^(1.0/3.0));              % Number of atoms that can fit on one side of cube
particleID = 1;                     % ID referring to particle in system
dx = side / n;                      % Separation space between particles (Angstrom / atom)

% For loop assigning positions to particles in perfect cube configuration.
% Later on we will use data to set initial poisitions.
for x = 1 : n
    for y = 1 : n
        for z = 1 : n
            if (particleID <= N)
                % Set positions of particle i as (x, y, z) in position
                % matrix.
                particle = Particle("H", N, particleID);                % Instantiate a new Hydorgen particle
                particle.Position = [(dx * x) (dx * y) (dx * z)];       % Update the position of this particle
                pos(particleID) = particle;                             % Add particle's information into array
                particleID = particleID + 1;
            end
        end
    end
end
posMatrix = pos;
end
