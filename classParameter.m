% Definition of Object: classParameter
% Dario Colombari - 24/07/2023

classdef classParameter
    properties (Constant)

        n  {mustBeInteger, mustBePositive} = 1;
       
    end
        
    properties (Access = public)

        electricalLoadMultiplier  {mustBeNonnegative}   = 0;

    end
        
end
