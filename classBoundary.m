% Definition of Object: classBoundary
% Dario Colombari - 23/08/2023

classdef classBoundary
    properties (Constant)

        n  {mustBeInteger, mustBePositive} = 5;
       
    end
        
    properties (Access = public)

        Portata_Biogas      {mustBeNonnegative}          = 0;
        Conc_CH4            {mustBeNonnegative}          = 0;
        P_Carico            {mustBeReal}                 = 0;
        Prezzo              {mustBeReal}                 = 0;
        P_fv                {mustBeReal}                 = 0;

        % Boundaries for 
        serviceLife         {mustBeNonnegative}          = 20.0;



    end
        
end
