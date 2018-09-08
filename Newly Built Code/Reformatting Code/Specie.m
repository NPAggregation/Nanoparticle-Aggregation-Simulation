%% Specie Class %%
% This data structure holds chemical and physical information on a
% particle, atom, element, compound and etc.

classdef Specie
    properties
       elementName;                         % Name of element
       elementWeight;                       % Molecular weight of specie
       elementCharge;                       % Charge of specie
    end
    methods
        % Construct for class that takes in element/compounds symbol
        function specie = Specie(name)
           if nargin == 1
               if isstring(name)
                   if strcmp(name, 'H')
                       specie.elementName = 'Hydrogen';
                       specie.elementWeight = 1;
                       specie.elementCharge = 0;
                   end
               else
                   error('Parameter must be a string.')
               end
           end
        end
    end
end