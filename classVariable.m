% Definition of Object: classVariable
% Dario Colombari - 24/07/2023

classdef classVariable
    properties (Constant)

        n  {mustBeInteger, mustBePositive} = 7;
       
    end
        
    properties (Access = public)

        PnomFV                  {mustBeNonnegative}                 = 0;
        Dimensione_Stoccaggio   {mustBeNonnegative}                 = 0;
        n_Elettrolizzatori      {mustBeInteger, mustBeNonnegative}  = 0;
        Numero_Stack            {mustBeInteger,mustBeNonnegative}   = 0;
        Coeff_importazione_el   {mustBeReal}                        = 0;
        Coeff_accensione_el     {mustBeReal}                        = 0;
        Coeff_accensione_FC     {mustBeReal}                        = 0;

    end
        
end
