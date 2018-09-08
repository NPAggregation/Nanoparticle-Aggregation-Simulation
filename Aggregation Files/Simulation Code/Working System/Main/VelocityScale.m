%% Velocity and Accelration Scaling %%
% Scale the velocity and acceleration of particles based on temperature
% correction factor for the system.

function particles = VelocityScale(N, Particle, Tcorr)
for i = 1:N
    sumVSq = sum(Particle(i).Velocity.^2);
    sumASq = sum(Particle(i).Acceleration.^2);
    facV = sqrt(Tcorr / sumVSq);
    facA = sqrt(Tcorr / sumASq);
    Particle(i).Velocity = Particle(i).Velocity * facV;
    Particle(i).Acceleration = Particle(i).Acceleration * facA;
end
particles = Particle;
end