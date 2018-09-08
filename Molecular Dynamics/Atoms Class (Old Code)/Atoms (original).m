%%Class for Atoms
classdef Atoms
    properties
        atom
        protonNum
    end
    methods
        function obj = Atoms(a)
            if nargin > 0
                a = lower(a);
                switch a
                    case 'carbon'
                        obj.atom = 'C';
                        obj.protonNum = 6;
                    case 'oxygen' 
                        obj.atom = 'O';
                        obj.protonNum = 8;
                    case 'nitrogen'
                        obj.atom = 'N';
                        obj.protonNum = 7;
                    case 'gold'
                        obj.atom = 'Au';
                        obj.protonNum = 79;
                    case 'silver' 
                        obj.atom = 'Ag';
                        obj.protonNum = 49;
                    case 'hydrogen'
                        obj.atom = 'H';
                        obj.protonNum = 1;
                    case 'aluminum'
                        obj.atom = 'Al';
                        obj.protonNum = 13;
                    case 'sodium'
                        obj.atom = 'Na';
                        obj.protonNum = 11;
                    case 'magnesium'
                        obj.atom = 'Mg';
                        obj.protonNum = 12;
                    case 'iron'
                        obj.atom = 'Fe';
                        obj.protonNum = 56;
                    otherwise
                        disp('none of those atoms are in the system')
                end 
            end 
        end 
    end
end
                
        