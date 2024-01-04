function ICECase = ICESimulation(b, electricalLoadMultiplier) %#codegen

%% Initialization of Output Variables

ICECase.totalElectricity = 0;
ICECase.totalICEElectricity = 0;
ICECase.totalPurchasedElectricity = 0;
ICECase.totalSoldElectricity = 0;
ICECase.totalSelfElectricity = 0;
ICECase.totalCost = 0;
ICECase.energyPrice = 0;

%% Time Interval

% Service Life
b.serviceLife = 20;

%% Electric Load

% Electrical Load Profile
P_Carico = b.P_Carico * electricalLoadMultiplier;

% Preliminary Calculation
LHV_Biogas = 9.94 * b.Conc_CH4/100; % 9.94 LHV of Methane only  in [ kWh/Nm3 ]
Pin_Biogas = LHV_Biogas .* b.Portata_Biogas; % [ kW ]

PNomICE = (max(Pin_Biogas) * 0.26814) ^ 1.069518717;
ICEEffNom =  0.26814 * PNomICE^0.065;

Sovrapprezzo_acquisto = 10;
Prezzo_acquisto = b.Prezzo + Sovrapprezzo_acquisto;

%% Costs 

% Purchasing cost
ICECost = (15648 * PNomICE^-0.536 + 830.23 * PNomICE^-0.268 + 0.17053 * PNomICE^-0.478) * PNomICE;
ICECostYear = ICECost / b.serviceLife;

ICEEleEfficiency = 0.26814*PNomICE^0.065;

Pin_Biogas_Adim = Pin_Biogas / max(Pin_Biogas);

ICEElePLEfficiency = -5.64656745091845 * Pin_Biogas_Adim.^5 + 18.8636413761601 * Pin_Biogas_Adim.^4 - 22.1300537765926* Pin_Biogas_Adim.^3 + 9.70809616295603 * Pin_Biogas_Adim.^2 - 0.0745682476757175 * Pin_Biogas_Adim + 0.279702413580693;

ICEEleEfficiencyCorrected = ICEEleEfficiency * ICEElePLEfficiency;


%% Yearly Cost/Revenues for 

ICEP = Pin_Biogas .* ICEEleEfficiencyCorrected;

soldPower = max(ICEP - P_Carico, 0);
purchasedPower = max(P_Carico - ICEP, 0);

RevenuesGrid = soldPower .* b.Prezzo / 1000.0;
CostGrid = purchasedPower .* Prezzo_acquisto / 1000.0;

ICECase.totalElectricity = sum(P_Carico) * b.serviceLife;
ICECase.totalICEElectricity = sum(ICEP) * b.serviceLife;
ICECase.totalPurchasedElectricity = sum(purchasedPower) * b.serviceLife;
ICECase.totalSoldElectricity = sum(soldPower)* b.serviceLife;
ICECase.totalSelfElectricity = ICECase.totalICEElectricity + ICECase.totalPurchasedElectricity - ICECase.totalSoldElectricity;
ICECase.totalCost = ICECost + (sum(CostGrid) - sum(RevenuesGrid)) * b.serviceLife;
ICECase.energyPrice = ICECase.totalCost * 1000 / ICECase.totalElectricity;

end
