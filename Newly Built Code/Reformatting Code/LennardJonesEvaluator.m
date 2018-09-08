%% Lennard Jones Evaluator %%
% Determine the Lennard Jones potential between each particle pair ij

function particle = LennardJonesEvaluator(Particle, N, sigma, eps)
sig6 = sigma^6; sig12 = sigma^12;
for i = 1:N
    for j = (i + 1):N
        distij = Particle(i).NeighborList(j);
        if (distij ~= 0)
            const = 4 * eps;
            collisionDiff = ((sig12 / distij^12) - (sig6 / distij^6));
            LJ = const * collisionDiff;
            Particle(i).LJ(Particle(j).id) = LJ;
            Particle(j).LJ(Particle(i).id) = LJ;
        end
    end
end
particle = Particle;
end