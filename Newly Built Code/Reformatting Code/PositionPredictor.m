%% New Position Predictor %%
% Generates new kinematic values for each particle

function particle = PositionPredictor(Particle, N, dtv)
% Compute new position, velocity, acceleration and more for each particle
for i = 1:N
   Particle(i).Position(1:3) = Particle(i).Position(1:3) + (Particle(i).Velocity(1:3) * dtv(1)) + (Particle(i).Acceleration(1:3) * dtv(2)) + (Particle(i).d3(1:3) * dtv(3)) + (Particle(i).d4(1:3) * dtv(4)) + (Particle(i).d5(1:3) * dtv(5));
   Particle(i).Velocity(1:3) = Particle(i).Velocity(1:3) + (Particle(i).Acceleration(1:3) * dtv(1)) + (Particle(i).d3(1:3) * dtv(2)) + (Particle(i).d4(1:3) * dtv(3)) + (Particle(i).d5(1:3) * dtv(4)); 
   Particle(i).Acceleration(1:3) = Particle(i).Acceleration(1:3) + (Particle(i).d3(1:3) * dtv(1)) + (Particle(i).d4(1:3) * dtv(2)) + (Particle(i).d5(1:3) * dtv(3));
   Particle(i).d3(1:3) = Particle(i).d3(1:3) + (Particle(i).d4(1:3) * dtv(1)) + (Particle(i).d5(1:3) * dtv(2));
   Particle(i).d4(1:3) = Particle(i).d4(1:3) + (Particle(i).d5(1:3) * dtv(1));
end
particle = Particle;
end