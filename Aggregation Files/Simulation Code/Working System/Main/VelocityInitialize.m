%% Initializing Velocity and Acceleration %%
% Initializes velocity and acceleration of all the particles in the system

function particles = VelocityInitialize(Particles, N, Tcorr)
% Assign Random Velocities and Acceleration (Between -1 and 1)
for i = 1:N
    Particles(i).Velocity(1) = (2 * rand) - 1;
    Particles(i).Velocity(2) = (2 * rand) - 1;
    Particles(i).Velocity(3) = (2 * rand) - 1;
    Particles(i).Acceleration(1) = (2 * rand) - 1;
    Particles(i).Acceleration(2) = (2 * rand) - 1;
    Particles(i).Acceleration(3) = (2 * rand) - 1;
end

sumV = zeros(1, 3);                                     % Sum of velocities in three axes
sumA = zeros(1, 3);                                     % Sum of acceleration in three axes;

% Initialize System with Zero Net Momentum
for i = 1:N
    sumV(1:3) = sumV(1:3) + Particles(i).Velocity(1:3);    % Get sum of velocities in the ith-direction
    sumA(1:3) = sumA(1:3) + Particles(i).Acceleration(1:3);
end
    
% Normally distributing velocity and acceleration of particles in the system
for i = 1:N
    Particles(i).Velocity(1:3) =  Particles(i).Velocity(1:3) - (sumV(1:3) / N);
    Particles(i).Velocity(1:3) =  Particles(i).Acceleration(1:3) - (sumA(1:3) / N);
end

% Velocity and Acceleration Scaling
for i = 1:N
    sumVSq = sum(Particles(i).Velocity.^2);
    sumASq = sum(Particles(i).Acceleration.^2);
    facV = sqrt(Tcorr / sumVSq);
    facA = sqrt(Tcorr / sumASq);
    Particles(i).Velocity = Particles(i).Velocity * facV;
    Particles(i).Acceleration = Particles(i).Acceleration * facA;
end

particles = Particles;
end