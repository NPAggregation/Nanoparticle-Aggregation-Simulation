%% Temperature Control %%
% Compute Temperature at a given moment in the system

function T = TempScaling(mass, N, N_adjust, kb, Particles)
total_velocity = 0;
for i = 1:N
    total_velocity = total_velocity + dot(Particles(i).Velocity, Particles(i).Velocity); 
end
T = (mass / (3 * N_adjust * kb)) * total_velocity;
end