%% Neighbour Class %%
% Stores the particle information of a neighbouring particles, as well as
% the distance from the particle of interest.

classdef Neighbour
    properties 
        Particle;                         % Particle data
        Distance;                         % Distance from particle of interest
    end
    methods
        % Constructor for neighbour
        function neighbour = Neighbour(particle, distance)
            if nargin == 2
               if isnumeric(distance)
                   neighbour.Distance = distance;
                   neighbour.Particle = particle;
               end
            end
        end
    end
end