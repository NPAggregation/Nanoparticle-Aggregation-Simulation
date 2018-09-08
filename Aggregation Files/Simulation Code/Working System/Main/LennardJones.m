%% Lennard Jones Calculator %%
% Calculate Lennard Jones Potential

function LJ = LennardJones(eps, sigma, rCut)
LJ = 4 * eps * ((sigma / rCut)^12 - (sigma / rCut)^6);  % Lennard Jones Potential (eV)
end