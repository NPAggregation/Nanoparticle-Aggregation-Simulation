%% Velocity Scale %%
% Function used to scale velocity of particles based on
% current and old temperature

function particles = Velocity_Scale(Particles, N, T0, T1)
for i = 1:N
    Particles(i).Velocity = Particles(i).Velocity * (sqrt(T0/T1));
end

particles = Particles;
end