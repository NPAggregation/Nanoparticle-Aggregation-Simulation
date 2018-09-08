%% Force Evaluator %%
% Determine the Lennard Jones force between each particle pair ij

function particle = ForceEvaluator(Particle, N, sigma, eps)
sig6 = sigma^6; sig12 = sigma^12;
for i = 1:N
    for j = (i + 1):N
        distij = Particle(i).NeighborList(j);
        if (distij ~= 0)
            const = (24 * eps) / distij;
            collisionDiff = ((2 * (sig12 / distij^12)) - (sig6 / distij^6));
            force = const * collisionDiff;
            Particle(i).Force(Particle(j).id) = force;
            Particle(j).Force(Particle(i).id) = -force;
        end
    end
end
particle = Particle;
end