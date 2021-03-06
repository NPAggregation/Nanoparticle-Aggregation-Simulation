%% New Position Predictor %%
% Generates new kinematic values for each particle
% Using Verlet Integration formulae's

function particle = PositionPredictor(Particle, N, dt, mass)
rmass = 1.0 / mass;
% Compute new position, velocity, acceleration and more for each particle
for i = 1:N
   Particle(i).Position(1:3) = Particle(i).Position(1:3) + (Particle(i).Velocity(1:3) * dt) + (0.5 * Particle(i).Acceleration(1:3) * dt * dt);
   Particle(i).Velocity(1:3) = Particle(i).Velocity(1:3) + (0.5 * (Particle(i).Force(1:3) * rmass + Particle(i).Acceleration(1:3)) * dt);
   Particle(i).Acceleration(1:3) = Particle(i).Force(1:3) * rmass * 10^15;
end
particle = Particle;
end