%% System Temperature Calculator %%
% Uses only the kinetic energy of the system to determine temperature.
% A modified form of the ideal gas constant, with adjusting parameters.
% Used as part of equilibration process.

function T = ComputeSystemTemperature(KE, Ek_Adjust, Na, R)
T = (Ek_Adjust * Na * KE) / R;
end