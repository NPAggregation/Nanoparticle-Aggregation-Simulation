%% Temperature Scaling Function %%
% Function is used to adjust the temperature of the system to allow
% equlibration. Bounds have been set on the value of the constant, as not
% to cause INFINITY error in particle motion parameters.

function velocity_scale = TemperatureScale(T, T0, bounds)
velocity_scale = sqrt(T/T0);
if (velocity_scale > bounds(1))
    velocity_scale = bounds(1);
elseif (velocity_scale < bounds(2))
    velocity_scale = bounds(2);
end
end