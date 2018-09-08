%% Repel Particles %%
% Used to avoid overlapping of particles, and create a repulsion effect.

function particles = RepelParticles(Particles, N)
responseFactor = 0.25;          % Randomly determined constant for repulsion effect
cutOffDist = 3;                 % Randomly determined cut-off distance
for i = 1:N
   for j = (i + 1) : N
      distance = Particles(i).NeighborList(j);
      % If distance threshold is met, repel particle in opposite direction
      if (distance <= cutOffDist)
         Particles(i).Acceleration = -(distance / responseFactor) * Particles(i).Acceleration;
      end
   end
end
particles = Particles;
end