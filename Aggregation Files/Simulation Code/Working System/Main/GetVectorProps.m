%% Get Vector Properties %%
% Used to get a particular property from particles data set.
% Property could be velocity, acceleration, etc.

function vector = GetVectorProps(particles, N, pos)
variable = zeros(1, N);
for i = 1:N
   variable(i) = particles(i).Position(pos); 
end
vector = variable;
end