%% Boundary Condition %%
% Checks that particles haven't left system. If they have, they're forced
% back into the simulation space.

function particles = BoundaryCondition(N, Particle, side)
for i = 1:N
   for j = 1:3
      if (Particle(i).Position(j) > side); Particle(i).Position(j) = Particle(i).Position(j) - side; end
      if (Particle(i).Position(j) < 0); Particle(i).Position(j) = Particle(i).Position(j) + side; end
   end
end
particles = Particle;
end