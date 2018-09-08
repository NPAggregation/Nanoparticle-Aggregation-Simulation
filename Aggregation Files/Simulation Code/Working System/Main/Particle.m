%% Particle Class %%
% The data structure holds the following information about a particle in
% the system. The list of neighbours the particle of interest has. The
% position of the particle in our grid as well as other mechanical properties.
% Then there's a class of species which holds the element/compound of the 
% specie, containing it's atomic weight, charge and other relevant data to 
% a specie.

classdef Particle
    properties
        SpecieData;                     % Holds information on specie
        NeighborList;                   % List of neighbouring particles
        Position;                       % Position of particle in system
        Velocity;                       % Velocity Matrix (Angstrom/s)
        Acceleration;                   % Acceleration Matrix (Angstrom/s^2)
        d3;                             % Third Derivative (Snap) (Angstrom/s^3)
        d4;                             % Fourth Derivative (Crackle) (Angstrom/s^4)
        d5;                             % Fifth Derivative (Crackle) (Angstrom/s^5)
        Force;                          % Force ((Grams Angstrom)/s^2)
        LJ;                             % Lennard Jones Potential (eV)
        Momentum;                       % Momentum ((Grams Angstrom)/s)
        id;                             % Particle Identifier Number 
    end
    methods
        % Constructor of particle data class
        % Takes in the symbol of element/compound and size of neighbourlist
        function particle = Particle(specieSymbol, arrSize, id)
           if nargin == 3
              % Initialize the following arrays with x, y , z
              particle.Position = zeros(1, 3);
              particle.Velocity = zeros(1, 3);
              particle.Acceleration = zeros(1, 3);
              particle.d3 = zeros(1, 3);
              particle.d4 = zeros(1, 3);
              particle.d5 = zeros(1, 3);
              particle.Momentum = zeros(1, 3);
              particle.id = id;
              if isnumeric(arrSize)
                  particle.NeighborList = zeros(1, arrSize);       % Instantiate neigboour list with number of particles in system as worst case
                  particle.Force = zeros(1, 3);
                  particle.LJ = zeros(1, arrSize);
              else
                  error('Second parameter must be numeric.')
              end
              if isstring(specieSymbol)
                  particle.SpecieData = Specie(specieSymbol);
              else
                  error('First parameter must be a string')
              end
           end
        end
    end
end