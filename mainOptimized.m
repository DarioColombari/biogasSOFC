%% Script Initialization

close all
clear 
clc

p = classParameter;
v = classVariable;
b = classBoundary;

alreadyCompiled = false;

%% Setting Simulation Boundaries

% Reading Hourly Photovoltaic Plant Data [ W/kWNom ]
b.P_fv = table2array(readtable("Dati FV.csv", 'Range','B2:B8761'));

% Reading Hourly Biogas Flow Rates [ Nm3/h ]
b.Portata_Biogas = table2array(readtable('Dati Biogas.csv','Range','B2:B8761'));

% Reading Hourly Biogas Mass Concentrations [ - ]
b.Conc_CH4 = table2array(readtable('Dati Biogas.csv','Range','C2:C8761'));

% Reading Hourly Electricity Price [ €/kWh ]
b.Prezzo = table2array(readtable('Dati Prezzi.xlsx','Range','C2:C8761'));

% Reading Electricity Load by Energy Community [ kW ]
b.P_Carico = table2array(readtable('Dati Carico Elettrico.csv','Range','B2:B8761'));

%% Setting Simulation Boundaries

% Nominal Photovoltaic Installed Power [ kW ]
% PnomICE = 800;
% A = 16.4 * PnomICE + 4008;
% nomFVMax = 0.643 * A / 36;
% v.PnomFV = linspace(0, nomFVMax, 4);
v.PnomFV = [ 1500:500:6000 ];

% Hydrogen Storage Size [ Nm3 ]
v.Dimensione_Stoccaggio = [ 1000:1000:5000 ];

% Number of Electrolyser in Series [ - ]
v.n_Elettrolizzatori =  [ 0 5 10 15 20 ];

% Number of Fuel Cells
v.Numero_Stack = [  350:50:650 ];

% Threshold Value for Purchasing Electricity [ - ]
v.Coeff_importazione_el = [ 0 0.05 0.1 0.15 0.2  ];

% Threshold Value for Turning On The Electrolyser [ - ]
v.Coeff_accensione_el = [ 0 0.05 0.1 0.15 0.2 ];

% Threshold Value for Turning On The Fuel Cell [ - ]
v.Coeff_accensione_FC = [  0.8 1.0 1.2 1.4 1.6 ];

%% Setting Sensitivity Analysis Parameter

p.electricalLoadMultiplier = [0.4 : 0.2 : 1.8];

%% Compiling Solution

constantValue = 1;

if not(alreadyCompiled)
    
    cfg = coder.config('mex');
    cfg.IntegrityChecks = false;

    codegen -config cfg PPSimulation -args { b, v, constantValue}

    codegen -config cfg ICESimulation -args { b, constantValue}

    codegen -config cfg XSimulation -args { b, constantValue }
end



%% Running Simulation in Parallel with Compiled Code for Speed

t1 = tic;
parfor i = 1:length(p.electricalLoadMultiplier)

    optimizedCase(i) = PPSimulation_mex(b, v, p.electricalLoadMultiplier(i));
    ICECase(i) = ICESimulation_mex(b, p.electricalLoadMultiplier(i));
    XCase(i) = XSimulation_mex(b, p.electricalLoadMultiplier(i));
    
    disp("Simulation with Load Multiplier " + p.electricalLoadMultiplier(i) + " complete.")

end
t1 = toc(t1);

%% Results post processing

% Plotting Optimal Photovoltaic Size
figure("Name", "Optimal Photovoltaic Size")
plot(p.electricalLoadMultiplier, [optimizedCase.PnomFV]);
title("Optimal Photovoltaic Size")
xlim([min(p.electricalLoadMultiplier), max(p.electricalLoadMultiplier)])
xlabel("Load Multiplier {\it [-]}")
ylim([min(v.PnomFV) - 500, max(v.PnomFV) + 500])
ylabel("Optimal Photovoltaic Size {\it [kW]}")
grid on
saveas(gcf,'.\images\OptFV.svg');

% Plotting Optimal Storage Size
figure("Name", "Optimal H_2 Storage Size")
plot(p.electricalLoadMultiplier, [optimizedCase.Dimensione_Stoccaggio]);
title("Optimal H_2 Storage Size")
xlim([min(p.electricalLoadMultiplier), max(p.electricalLoadMultiplier)])
xlabel("Load Multiplier {\it [-]}")
ylim([min(v.Dimensione_Stoccaggio) - 500, max(v.Dimensione_Stoccaggio) + 500])
ylabel("Optimal H_2 Storage Size {\it [Nm^3]}")
grid on
saveas(gcf,'.\images\OptStorage.svg');

% Plotting Optimal Number of Electrolysers
figure("Name", "Optimal Number of Electrolysers")
plot(p.electricalLoadMultiplier, [optimizedCase.n_Elettrolizzatori]);
title("Optimal Number of Electrolysers")
xlim([min(p.electricalLoadMultiplier), max(p.electricalLoadMultiplier)])
xlabel("Load Multiplier {\it [-]}")
ylabel("Optimal Number of Electrolysers {\it [-]}")
grid on
saveas(gcf,'.\images\OptNElectrolysers.svg');

% Plotting Optimal Number of Fuel Cells in Stack
figure("Name", "Optimal Number of Fuel Cells in Stack")
plot(p.electricalLoadMultiplier, [optimizedCase.Numero_Stack]);
title("Optimal Number of Fuel Cells in Stack")
xlim([min(p.electricalLoadMultiplier), max(p.electricalLoadMultiplier)])
xlabel("Load Multiplier {\it [-]}")
ylabel("Optimal Number of Cells {\it [-]}")
grid on
saveas(gcf,'.\images\OptNFuelCell.svg');

% Plotting Optimal Electrolyser Energy Import Coefficient
figure("Name", "Optimal Electrolyser Energy Import Coefficient")
plot(p.electricalLoadMultiplier, [optimizedCase.Coeff_importazione_el]);
title("Optimal Electrolyser Energy Import Coefficient")
xlim([min(p.electricalLoadMultiplier), max(p.electricalLoadMultiplier)])
xlabel("Load Multiplier {\it [-]}")
ylabel("Optimal Electrolyser Import Factor {\it [-]}")
grid on
saveas(gcf,'.\images\OptImportCoeff.svg');

% Plotting Optimal Fuel Cell Energy Import Coefficient
figure("Name", "Optimal Fuel Cell Power On Coefficient")
plot(p.electricalLoadMultiplier, [optimizedCase.Coeff_accensione_FC]);
title("Optimal Fuel Cell Power On Coefficient")
xlim([min(p.electricalLoadMultiplier), max(p.electricalLoadMultiplier)])
xlabel("Load Multiplier {\it [-]}")
ylabel("Optimal Fuel Cell Power On Coefficient {\it [-]}")
grid on
saveas(gcf,'.\images\OptFCCoeff.svg');

% Plotting Optimal Fuel Cell Energy Import Coefficient
figure("Name", "Optimal Electrolyser Power On Coefficient")
fig=plot(p.electricalLoadMultiplier, [optimizedCase.Coeff_accensione_el]);
title("Optimal Electrolyser Power On Coefficient")
xlim([min(p.electricalLoadMultiplier), max(p.electricalLoadMultiplier)])
xlabel("Load Multiplier {\it [-]}")
ylabel("Optimal Electrolyser Power On Coefficient {\it [-]}")
grid on
saveas(fig,'.\images\OptEleCoeff.svg');

% Plotting Optimal Fuel Cell Energy Import Coefficient
energySummary = [-[optimizedCase.P_Carico_Tot]; [optimizedCase.P_FC_car_Tot]; [optimizedCase.P_FV_car_Tot]; -[optimizedCase.P_FC_rete_Tot]; -[optimizedCase.P_FV_rete_Tot]; [optimizedCase.P_importata_el_Tot]]';
figure("Name", "Energy Summary")
bar([-[optimizedCase.P_Carico_Tot]; [optimizedCase.P_FC_car_Tot]; [optimizedCase.P_FV_car_Tot]; -[optimizedCase.P_FC_rete_Tot]; -[optimizedCase.P_FV_rete_Tot]; [optimizedCase.P_importata_el_Tot]]', 'stacked');
title("Optimal Electrolyser Power On Coefficient")
xlabel("Load Multiplier {\it [-]}")
ylabel("Optimal Electrolyser Power On Coefficient {\it [-]}")
grid on
saveas(gcf,'.\images\EnergySummary.svg');

% Comparison of Total Costs
figure("Name", "Comparison of Solution")
bar([XCase.totalCost; ICECase.totalCost; optimizedCase.costoTotale ]', 'grouped');
title("Comparison of Solutions")
xlabel("Configuration {\it [-]}")
ylabel("Total Cost for Electricity {\it [€]}")
grid on
saveas(gcf,'.\images\costComparison.svg');

% Comparison of Total Specific Costs
figure("Name", "Comparison of Specific Cost for Solution")
bar([XCase.energyPrice; ICECase.energyPrice;  optimizedCase.P_Carico_Tot ./ optimizedCase.costoTotale]', 'grouped');
title("Comparison of Specific Cost for Different Solutions")
xlabel("Configuration {\it [-]}")
ylabel("Specific Total Cost for Electricity {\it [€/MWh]}")
grid on
saveas(gcf,'.\images\SpecificCostComparison.svg');

save("simulationResults.mat", "b", "p", "v", "optimizedCase", "XCase", "ICECase");

disp("Simulation Complete!")
