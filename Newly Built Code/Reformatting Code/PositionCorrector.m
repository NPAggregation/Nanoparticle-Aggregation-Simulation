%% Positions Corrector %%
% Function is to correct/adjust kinematic properties of each particle.
% Correction factor is determined based on the difference between the
% particles calculated acceleration and expected acceleration.

function particle = PositionCorrector(Particles, N, mass, alpha)
for i = 1:N
   errVec(1:3) = (sum(Particles(i).Force) / mass) - Particles(i).Acceleration(1:3);
   Particles(i).Position(1:3) = Particles(i).Position(1:3) + (errVec(1:3) * alpha(1));
   Particles(i).Velocity(1:3) = Particles(i).Velocity(1:3) + (errVec(1:3) * alpha(2));
   Particles(i).Acceleration(1:3) = Particles(i).Acceleration(1:3) + (errVec(1:3) * alpha(3));
   Particles(i).d3(1:3) = Particles(i).d3(1:3) + (errVec(1:3) * alpha(4));
   Particles(i).d4(1:3) = Particles(i).d4(1:3) + (errVec(1:3) * alpha(5));
   Particles(i).d5(1:3) = Particles(i).d5(1:3) + (errVec(1:3) * alpha(6));
end
particle = Particles;
end