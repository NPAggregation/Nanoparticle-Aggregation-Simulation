%% Total Velocity Calculator %%
% Calculator the total amount of velocity in the system.

function v = TotalVelocity(N, Particles)
v = 0;
for i = 1:N
   v = v + sum(Particles(i).Velocity); 
end
end