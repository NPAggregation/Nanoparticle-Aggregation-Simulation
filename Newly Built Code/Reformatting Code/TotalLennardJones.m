%% Total Lennard Jones Calculator %%
% Calculate the total Lennard Jones Potential in the system.

function U = TotalLennardJones(N, Particles)
U = 0;
for i = 1:(N - 1)
   for j = (i + 1):N
      U = U + Particles(i).LJ(j); 
   end
end
end