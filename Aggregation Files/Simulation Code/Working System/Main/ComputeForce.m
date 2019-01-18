%% Force Computing Function %%
% Determine the force experienced by each particle

function particle = ComputeForce(Particle, N, velocity_scale)
for i = 1:N
    for j = (i+1):N
        if (i ~= j)
            distance(1:3) = Particle(i).Position(1:3) - Particle(j).Position(1:3);
            totalDistance = sqrt(sum(distance.^2));
            truncDistance = min(totalDistance, pi / 2.0);
            
            Particle(i).Force(1:3) = (Particle(i).Force(1:3) - distance(1:3).*sin(2.0 * truncDistance) / totalDistance) / velocity_scale;
            Particle(j).Force(1:3) = (Particle(j).Force(1:3) - distance(1:3).*sin(2.0 * truncDistance) / totalDistance) / velocity_scale;
        end
    end
end
particle = Particle;
end