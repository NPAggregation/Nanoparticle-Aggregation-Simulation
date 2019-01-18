%% Initialize Velocity %%
% Function used to intialize velocity of all particles in system
% using Maxwell Boltzmann Distribution

function particle = InitializeVelocity(Particle, N, T)
% Assign random initial velocity
for i = 1:N
    r = -1 + rand(1,3) * (1 - -1);
    
    Particle(i).Velocity = r * sqrt(T);
end

% Minimize total linear momentum
for i = 1:N
    total_velocity = 0;
    for j = 1:N
        total_velocity = total_velocity + Particle(j).Velocity;
    end
    
    Particle(i).Velocity = Particle(i).Velocity - (total_velocity / N);
end

particle = Particle;
end