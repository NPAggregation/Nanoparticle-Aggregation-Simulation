%% Distance Computing Method %%
% Computes distance between each particle in the system
% Allows us to set the particle distance of (i,j) and (j,i)

function particle = ComputeDistance(Particle, N, rNbr2)
for i = 1:N
    for j = (i+1):N
       if (i ~= j)
            % Compute distance between ith and jth particle, where i doesnt equal j
            distance(1:3) = (Particle(i).Position(1:3) - Particle(j).Position(1:3)).^2;
            distance = sqrt(sum(distance));
            % If so update distance of both i and j particles
            if (distance <= rNbr2)
                neighborIDi = Particle(j).id;
                neighborIDj = Particle(i).id;
                Particle(i).NeighborList(neighborIDi) = distance;
                Particle(j).NeighborList(neighborIDj) = distance;
            end 
       end      
    end
end
particle = Particle;
end