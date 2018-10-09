%% Total Velocity Calculator %%
% Calculator the total amount of velocity in the system.

function v = TotalVelocity(N, Particles)
v = 0;
for i = 1:N
   v = v + sqrt(sum(Particles(i).Velocity(1:3).^2)); 
end
end