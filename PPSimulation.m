function opt = PPSimulation(b, v, electricalLoadMultiplier) %#codegen

%% Initialization of optimization results

opt.PnomFV = 0;
opt.n_Elettrolizzatori = 0;
opt.Dimensione_Stoccaggio = 0;
opt.Numero_Stack = 0;
opt.Coeff_importazione_el = 0;
opt.Coeff_accensione_el = 0;
opt.Coeff_accensione_FC = 0;
opt.P_Carico_Tot = 0;
opt.P_FC_car_Tot = 0;
opt.P_FV_car_Tot = 0;
opt.P_FC_rete_Tot = 0;
opt.P_FV_rete_Tot = 0;
opt.P_importata_el_Tot = 0;
opt.costoTotale = 0;

% Yearly Simulation Resolution
DeltaT = 1; % [ h ]
endTime = 8760 / DeltaT; % Total Number of Timesteps

%% Optimization Variables

% Number of Total Simulate Cases and Simulation Progress
NTotalCases = length(v.PnomFV)*length(v.Dimensione_Stoccaggio)*length(v.n_Elettrolizzatori)*length(v.Numero_Stack)*length(v.Coeff_accensione_FC)*length(v.Coeff_importazione_el)*length(v.Coeff_accensione_el);
P_CaricoBase = b.P_Carico;

caseCounter = 0;
P_Carico = P_CaricoBase * electricalLoadMultiplier;

%% Preliminary Calculations


% Costs of Photovoltaic Plant
CostoINV_FV = 684 * v.PnomFV; % 684 da scenari costi futuri IRENA
VitaUtileFV = 20;
CostoINV_FV_annuo = CostoINV_FV / VitaUtileFV;
CostoOEM_FV = 0.02 * CostoINV_FV; % costo annuo

% Operations of Anaerobic Digester
LHV_Biogas = 9.94 * b.Conc_CH4/100; % 9.94 LHV of Methane only  in [ kWh/Nm3 ]
Pin_Biogas = LHV_Biogas .* b.Portata_Biogas; % [ kW ]

% Operations of Reformer
PortataH2_Biogas = 1.4902 * b.Portata_Biogas; % qui andranno inseriti i risultati del modello reformer 0.8367

% Costs of Reformer
% No information available now.

% Costs of Anaerobic Digester
%   Useful for a correct estimation of LCOE of the plant. Useless for
%   system optimization and short-term operations as the digester has no
%   avoidable costs.

% Operations of Electrolyser
Pin_singolo_Elettrolizzatore = 2.7; % Potenza nominale in kW in ingresso al singolo elettrolizzatore
Pin_nominale_Elettrolizzatore = Pin_singolo_Elettrolizzatore * v.n_Elettrolizzatori;
Portata_H2_El_Nominale = 0.5 * v.n_Elettrolizzatori; % Portata nominale elettrolizzatore in Nm3/h

% Costs of Electrolyser
CostoINV_elettrolizzatore = 200 * Pin_nominale_Elettrolizzatore; % 667 da scenari costi futuri DOE
VitaUtileElettrolizzatore = 10;
CostoINV_elettrolizzatore_annuo = CostoINV_elettrolizzatore / VitaUtileElettrolizzatore;
CostoOEM_elettrolizzatore = 0.03 * CostoINV_elettrolizzatore; % costo annuo

% Operations of Storage

% Cost of Hydrogen Storage
CostoINV_stoccaggio = v.Dimensione_Stoccaggio * 40; % 40 €/Nm3 da dati scambiati con Luigi
VitaUtileStoccaggio = 20;
CostoINV_stoccaggio_annuo = CostoINV_stoccaggio / VitaUtileStoccaggio;
CostoOEM_stoccaggio = 0.01 * CostoINV_stoccaggio; % costo annuo

% Operations of Fuel Cell
PotenzaNominale_Stack = 1.77; % kW di potenza nominale stack FC
LitriNom_Stack = 57.728; % L/min di H2 del singolo stack FC
PotenzaNominale_FC = PotenzaNominale_Stack * v.Numero_Stack; % Potenza massima FC
LitriNom_FC = LitriNom_Stack * v.Numero_Stack; %Portata massima inviabile alla FC
Temperatura_FC = 1023; %T stack in K
Pressione_FC = 1.2; %p stack in bar
thermal_eff_FC = 0.31; %thermal efficiency of the FC
PortataNominale_FC = LitriNom_Stack * v.Numero_Stack / 1000 * 60 * Pressione_FC / 1.01325 * 273.15 / Temperatura_FC; %in Nm3/h

% Costs of Fuel Cell
CostoINV_FC = 2500 * PotenzaNominale_FC; % 3340 Target Cosmos usato da Santarelli & Co
VitaUtileFC = 10;
CostoINV_FC_annuo = CostoINV_FC / VitaUtileFC;
CostoOEM_FC = 0.01 * CostoINV_FC; % costo annuo

% Costs of Electricity Market
Prezzo_medio = mean (b.Prezzo);
Sovrapprezzo_acquisto = 10;
Prezzo_acquisto = b.Prezzo + Sovrapprezzo_acquisto;
Prezzo_acquisto_medio = mean (Prezzo_acquisto);

endTime = endTime + 1;
%% Preallocation of Variables

Pin_El0 = zeros(endTime,1);
Pin_El = zeros(endTime,1);
PortataH2_El = zeros(endTime,1);
P_impiantoFV = zeros(endTime,1);
P_FV_rete = zeros(endTime,1);
P_FV_car = zeros(endTime,1);
Livello_stoccaggio = zeros(endTime,1);
Livello_stoccaggio0 = zeros(endTime, 1);
portata_H2_FC = zeros(endTime, 1);
portata_H2_FC0 = zeros(endTime, 1);
litri_H2_FC = zeros(endTime, 1);
litri_H2_FC0 = zeros(endTime, 1);
P_FC = zeros(endTime, 1);
P_FC_car = zeros(endTime, 1);
P_FC_car0 = zeros(endTime, 1);
P_FC_rete = zeros(endTime, 1);
P_importata_el = zeros(endTime, 1);
P_importata = zeros(endTime, 1);
P_in_FC = zeros(endTime, 1);
P_th_FC = zeros(endTime, 1);
P_th_reformer = zeros(endTime, 1);
P_el_reformer = zeros(endTime, 1);
Costo_P_importata = zeros(endTime, 1);
kk=1;
opt_PnomFV = 0;
opt_n_Elettrolizzatori = 0;
opt_Dimensione_Stoccaggio = 0;
opt_Numero_Stack = 0;
opt_Coeff_importazione_el = 0;
opt_Coeff_accensione_el = 0;
opt_Coeff_accensione_FC = 0;
minimo = 1/eps;
Livello_inizio_anno = 0;
Costo_P_importata_annuo = zeros(NTotalCases, 20);
Costo_totale = zeros(NTotalCases, 1);
Costo_P_importata_totale = zeros(NTotalCases, 1);
endTime = endTime - 1;

%% Yearly Simulation for the Different Configurations

% Optimization Parameters 
for a = 1:length (v.PnomFV)

  for h = 1:length (v.n_Elettrolizzatori)
     
     for c = 1:length (v.Dimensione_Stoccaggio)
         
         for d = 1:length (v.Numero_Stack)

             for e = 1:length (v.Coeff_importazione_el)

                 for f = 1:length (v.Coeff_accensione_el)

                     for g = 1:length (v.Coeff_accensione_FC)
                       
                        P_Carico_Tot = 0;
                        P_FC_car_Tot = 0;
                        P_FV_car_Tot = 0;
                        P_FC_rete_Tot = 0;
                        P_FV_rete_Tot = 0;
                        P_importata_el_Tot = 0;

                        % Simulation for Whole Plant Life 
                        for j=1:20

                        Livello_stoccaggio(1) = 0; %Livello_inizio_anno;
                        Livello_stoccaggio0(1) = 0; %Livello_inizio_anno;

                        % Yearly Simulation
                        for i=1:DeltaT:8760
                     
                         P_th_reformer (i) = 52.4 / 49.21 * b.Portata_Biogas(i); %kW termici consumati dal reformer per 49.21 Nm3/h di Biogas; inserire dati da modello reformer
                         P_el_reformer (i) = 7.5 /49.21 * b.Portata_Biogas(i); %kW elettrici consumati dal reformer per 49.21 Nm3/h di Biogas; inserire dati da modello reformer
                         P_impiantoFV (i) = b.P_fv(i) * v.PnomFV(a) /1000;

                       
% Logiche di funzionamento impianto: elettrolizzatore  
if Prezzo_acquisto(i) <= v.Coeff_importazione_el(e) * Prezzo_acquisto_medio %valore decisionale
    
     % Se il prezzo di acquisto è basso e il fotovoltaico è superiore ai carichi
    if P_impiantoFV(i) >= P_Carico(i) + P_el_reformer (i)
        P_FV_rete(i) = P_impiantoFV(i) - P_Carico(i) - P_el_reformer (i); % Nota: non dovrebbe esserci anche la potenza dell'elettrolizzatore?
    else % Se il prezzo è basso e il fotovoltaico è inferiore ai carichi
        P_FV_rete(i) = 0;
    end

    % Se prezzo di acquisto è basso, ma lo storage è pieno
    if Livello_stoccaggio(i) > 0.95 * v.Dimensione_Stoccaggio (c)
        Pin_El(i) = 0;
        PortataH2_El(i) = 0;
    else % Se prezzo di acquisto è basso e lo storage non è pieno
        Pin_El(i) = Pin_nominale_Elettrolizzatore (h);
        PortataH2_El(i) = Portata_H2_El_Nominale(h);
        if P_FV_rete(i) >= Pin_El(i)
            P_FV_rete(i) = P_FV_rete(i) - Pin_El(i);
            P_importata_el(i) = 0;
        else
            P_importata_el(i) = Pin_El(i) - P_FV_rete(i);
            P_FV_rete(i) = 0;
        end
    end

    % Se il prezzo di acquisto è alto e il fotovoltaico è superiore ai carichi
elseif P_impiantoFV(i) >= P_Carico(i) + P_el_reformer (i)
    
    % Se il prezzo di vendita è alto
    if b.Prezzo(i) >= v.Coeff_accensione_el(f) * Prezzo_medio %valore decisionale dovuto ad efficienza media elettrolizzatore+FC
        Pin_El0(i) = 0;
        PortataH2_El(i) = 0;
    else % Se il prezzo di vendita è basso
        Pin_El0(i) = P_impiantoFV(i) - P_Carico(i) - P_el_reformer (i);
    end
    
    % Se lo storage è pieno e il fotovoltaico non produce
    if Pin_El0(i) <= 0.2 * Pin_nominale_Elettrolizzatore (h) || Livello_stoccaggio(i) > 0.95 * v.Dimensione_Stoccaggio (c)
        Pin_El(i) = 0;
        PortataH2_El(i) = 0;
        P_FV_rete(i) = P_impiantoFV(i) - P_Carico(i) - P_el_reformer (i);
    elseif Pin_El0(i) <= Pin_nominale_Elettrolizzatore (h) && Pin_El0(i) > 0.2 * Pin_nominale_Elettrolizzatore (h) % Se il fotovoltaico produce sotto il nominale
        Pin_El(i) = Pin_El0(i);
        PortataH2_El(i) = (0.000180 *  Pin_El(i) * 1000 / v.n_Elettrolizzatori(h) + 0.013514) * v.n_Elettrolizzatori (h); %inserire curva elettrolizzatore
        P_FV_rete(i) = 0;
    else % Se il fotovoltaico produce più del nominale dell'elettrolizzatore
        Pin_El(i) = Pin_nominale_Elettrolizzatore (h); 
        PortataH2_El(i) = Portata_H2_El_Nominale(h);
        P_FV_rete(i) = P_impiantoFV(i) - P_Carico(i) - P_el_reformer (i) - Pin_El(i);
    end

    % Se il prezzo di acquisto è alto e il fotovoltaico è inferiore ai carichi 
else 
    Pin_El(i) = 0;
    PortataH2_El(i) = 0;
    P_FV_rete(i) = 0;
end

P_FV_car(i) = P_impiantoFV(i) - Pin_El(i) - P_FV_rete(i) + P_importata_el (i);

%Logiche di funzionamento impianto: FC

% Se lo storage è pieno
if Livello_stoccaggio(i) >= 0.95 * v.Dimensione_Stoccaggio(c)
    portata_H2_FC(i) = PortataNominale_FC(d);
    P_FC(i) = PotenzaNominale_FC(d);
    
    % Se la potenza di fotovoltaico e fuel cell è superiore ai carichi
    if P_FC(i) + P_FV_car(i) >= P_Carico(i) + P_el_reformer (i)
        P_FC_car(i) = P_Carico(i) + P_el_reformer (i) - P_FV_car(i);
        P_FC_rete(i) = P_FC(i)- P_FC_car(i);
    else % Se la potenza della fuel cell è inferiore ai carichi
        P_FC_rete(i) = 0;
        P_FC_car(i) = P_FC(i);
    end
    
 % Se il prezzo di vendita è elevato e lo storage non è pieno    
elseif b.Prezzo(i) >= v.Coeff_accensione_FC(g) * Prezzo_medio && Livello_stoccaggio(i) >= PortataNominale_FC(d) %valore decisionale
    portata_H2_FC(i) = PortataNominale_FC(d);
    P_FC(i) = PotenzaNominale_FC(d);

    % Se la potenza di fotovoltaico e fuel cell sono superiori al carico
    if P_FC(i) + P_FV_car(i) >= P_Carico(i) + P_el_reformer (i)
        P_FC_car(i) = P_Carico(i) + P_el_reformer (i) - P_FV_car(i);
        P_FC_rete(i) = P_FC(i) - P_FC_car(i);
    else % Se la potenza di fotovoltaico e fuel cell sono inferiori al carico
        P_FC_rete(i) = 0;
        P_FC_car(i) = P_FC(i);
    end

% Se il prezzo di vendita è elevato e lo storage è compreso tra 30-100 %
% della portata nominale della fuel cell
elseif b.Prezzo(i) >= v.Coeff_accensione_FC(g) * Prezzo_medio && Livello_stoccaggio(i) < PortataNominale_FC(d) && Livello_stoccaggio(i) >= 0.3 * PortataNominale_FC(d)
    portata_H2_FC(i) = Livello_stoccaggio(i);
    litri_H2_FC (i) = (portata_H2_FC(i) * 1.01325 / Pressione_FC * Temperatura_FC / 273.15) * 1000 / 60 / v.Numero_Stack(d); %portata in L/min per singolo stack in ingresso alla FC
    P_FC(i) = (-0.911778 * (litri_H2_FC (i) / LitriNom_Stack) ^ 2 + 2.400667 * (litri_H2_FC (i) / LitriNom_Stack) - 0.495745) * PotenzaNominale_FC(d); %inserire curva modello FC
    
    % Se la potenza di fotovoltaico e fuel cell sono inferiori al carico
    if P_FC(i) + P_FV_car(i) >= P_Carico(i) + P_el_reformer (i)
        P_FC_car(i) = P_Carico(i) + P_el_reformer (i) - P_FV_car(i);
        P_FC_rete(i) = P_FC(i) - P_FC_car(i);
    else % Se la potenza di fotovoltaico e fuel cell sono inferiori al carico 
        P_FC_rete(i) = 0;
        P_FC_car(i) = P_FC(i);
    end
else % Se lo storage non è pieno e il prezzo di vendita è basso
    P_FC_rete(i) = 0;
    P_FC_car0(i) = min (P_Carico(i) + P_el_reformer (i) - P_FV_car(i), PotenzaNominale_FC(d));
    
    %Se la potenza della fuel cell è sotto il 30%
    if P_FC_car0(i) < 0.3 * PotenzaNominale_FC(d)
        portata_H2_FC(i) = 0;
            P_FC_car(i) = 0;
            P_FC(i) = 0;
    else  %Se la potenza della fuel cell è sopra il 30%
        litri_H2_FC0(i) = (0.600489 * (P_FC_car0(i) / PotenzaNominale_FC(d)) ^2 + 0.080074 * (P_FC_car0(i) / PotenzaNominale_FC(d)) + 0.318264) * LitriNom_Stack; %inserire curva modello FC; in L/min
        portata_H2_FC0(i) = (litri_H2_FC0(i) * Pressione_FC / 1.01325 * 273.15 / Temperatura_FC) * 60 / 1000 * v.Numero_Stack(d); %in Nm3/h
        
        % Se la potenza della fuel cell era superiore al 30%
        if portata_H2_FC0(i) >= 0.3 * PortataNominale_FC (d)
            %Se il livello di stoccaggio era superiore alla portata  
            if Livello_stoccaggio(i) >= portata_H2_FC0(i)
               P_FC_car(i) = P_FC_car0(i);
               P_FC(i) = P_FC_car(i);
               portata_H2_FC(i) = portata_H2_FC0(i);

            % Se il livello di stoccaggio è maggiore al 30% della portata
            % della fuel cell
            elseif Livello_stoccaggio(i) >= 0.3 * PortataNominale_FC(d) 
                portata_H2_FC(i) = Livello_stoccaggio(i);
                litri_H2_FC (i) = (portata_H2_FC(i) * 1.01325 / Pressione_FC * Temperatura_FC / 273.15) * 1000 / 60 / v.Numero_Stack(d);
                P_FC_car(i) = (-0.911778 * (litri_H2_FC (i) / LitriNom_Stack) ^ 2 + 2.400667 * (litri_H2_FC (i) / LitriNom_Stack) - 0.495745) * PotenzaNominale_FC(d); %inserire curva modello FC
                P_FC(i) = P_FC_car(i);
            else % Se il livello stoccaggio è minore del 30% della portata 
                portata_H2_FC(i) = 0;
                P_FC_car(i) = 0;
                P_FC(i) = 0;
            end
        else % Se la portata della fuel cell era inferiore al 30%
            portata_H2_FC(i) = 0;
            P_FC_car(i) = 0;
            P_FC(i) = 0;
        end
    end
end

% Se la potenza della fuel cell e del fotovoltaico non sostengono i carichi
if  P_FC_car(i) + P_FV_car(i) >= P_Carico(i) + P_el_reformer (i)
    P_importata(i) = 0 +  P_importata_el (i) ;
else
    P_importata(i) = P_Carico(i) + P_el_reformer(i) - P_FC_car(i) - P_FV_car(i) +  P_importata_el (i);
end

P_in_FC (i) = portata_H2_FC(i) * 3; %P in FC in KW
P_th_FC (i) = P_in_FC (i) * thermal_eff_FC; %P th FC
Costo_P_importata (i) = P_importata (i) / 1000 * b.Prezzo(i); % costo annuo dell'elettricità importata

% Update Storage Level
Livello_stoccaggio0(i+1) = Livello_stoccaggio0(i) - portata_H2_FC(i)  + PortataH2_Biogas(i) + PortataH2_El(i);
if Livello_stoccaggio0(i+1) >= v.Dimensione_Stoccaggio(c)
    Livello_stoccaggio(i+1) = v.Dimensione_Stoccaggio(c);
else
    Livello_stoccaggio(i+1) = Livello_stoccaggio0(i+1);
end

%Energy Summary 
P_Carico_Tot = P_Carico_Tot + P_Carico(i);
P_FC_car_Tot = P_FC_car_Tot + P_FC_car(i);
P_FV_car_Tot = P_FV_car_Tot + P_FV_car(i);
P_FC_rete_Tot = P_FC_rete_Tot + P_FC_rete(i);
P_FV_rete_Tot = P_FV_rete_Tot + P_FV_rete(i);
P_importata_el_Tot = P_importata_el_Tot + P_importata(i);

end

% End of Year Update
if endTime + 1 <= 8760
    Livello_inizio_anno = Livello_stoccaggio(endTime+1); 
end
Costo_P_importata_annuo(kk,j) = sum(Costo_P_importata (:));
end



% Costs and Profits
Costo_P_importata_totale(kk) = sum(Costo_P_importata_annuo(kk,:));
Costo_totale(kk) = 1 * Costo_P_importata_totale(kk) + 20 * (CostoINV_FV_annuo (a) + CostoOEM_FV (a) + CostoINV_elettrolizzatore_annuo (h) + CostoOEM_elettrolizzatore (h) + CostoINV_stoccaggio_annuo (c) + CostoOEM_stoccaggio (c) + CostoINV_FC_annuo (d) + CostoOEM_FC (d));



             if Costo_totale(kk) < minimo
                % Save configuration if it is the optimum
                minimo = Costo_totale(kk);
                opt.PnomFV = v.PnomFV(a);
                opt.n_Elettrolizzatori = v.n_Elettrolizzatori(h);
                opt.Dimensione_Stoccaggio = v.Dimensione_Stoccaggio(c);
                opt.Numero_Stack = v.Numero_Stack(d);
                opt.Coeff_importazione_el = v.Coeff_importazione_el(e);
                opt.Coeff_accensione_el = v.Coeff_accensione_el(f);
                opt.Coeff_accensione_FC = v.Coeff_accensione_FC(g);
                opt.P_Carico_Tot = P_Carico_Tot;
                opt.P_FC_car_Tot = P_FC_car_Tot;
                opt.P_FV_car_Tot = P_FV_car_Tot;
                opt.P_FC_rete_Tot = P_FC_rete_Tot;
                opt.P_FV_rete_Tot = P_FV_rete_Tot;
                opt.P_importata_el_Tot = P_importata_el_Tot;
                opt.costoTotale = Costo_totale(kk);
             end
             kk = kk+1;


            end
        end
    end
end
             end
         end
     end




end
