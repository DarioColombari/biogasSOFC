function XCase = XSimulation(b, electricalLoadMultiplier) %#codegen

%% Initialization of Output Variables

XCase.totalElectricity = 0;
XCase.totalCost = 0;
XCase.energyPrice = 0;

%% Electric Load

% Electrical Load Profile
P_Carico = b.P_Carico * electricalLoadMultiplier;

%% Yearly Cost of Electricity from Grid

% Service Life (of equivalent power plant)
b.serviceLife = 20;

% Costs of Electricity Market
purchasingCost = b.Prezzo .* P_Carico / 1000;

%% Total Cost and Electricity, Average Electricity Price

% Total Electricity
XCase.totalElectricity = sum(P_Carico) * b.serviceLife;

% Total Cost
XCase.totalCost = sum(purchasingCost) * b.serviceLife;

% Average Electricity Price
XCase.energyPrice = XCase.totalCost * 1000 / XCase.totalElectricity;

end
