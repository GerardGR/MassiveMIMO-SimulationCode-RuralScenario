function [PL_CIH] = pathloss_CIH_v1(d, fc, h_BS)
%%
%% distance_  in m
%% fc  in GHz
%% h_BS: in m
%% derived from ICC_Rma, Equations 15 & 16, pp 6
%% shadowing fading std  LoS: 1.7 dB
%%                       NLoS: 6.7 dB
%%  
    X_sigma_LoS = 0;        %% Shadowing fading LoS
    X_sigma_NLoS = 0;       %% Shadowing fading NLoS

    %% Pathloss for LoS transmission
    PL_LoS = 32.4+20*log10(fc)+23.1*(1-0.03*((h_BS-35)/35))*log10(d)+ X_sigma_LoS;
    
    %% Pathloss for NLoS transmission 
    PL_NLoS = 32.4+20*log10(fc)+30.7*(1-0.049*((h_BS-35)/35))*log10(d)+ X_sigma_NLoS;
    
    %% Patloss function
    PL_CIH = probability_LoS(d)*PL_LoS + (1-probability_LoS(d))*PL_NLoS;