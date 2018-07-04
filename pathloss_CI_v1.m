function [PL_CI] = pathloss_CI_v1(d, fc)
%%
%% distance_  in m
%% fc  in GHz
%% derived from ICC_Rma, Equations 13 & 14, pp 6
%% shadowing fading std  LoS: 1.7 dB
%%                       NLoS: 6.7 dB
%%
    X_sigma_LoS = 0;        %% Shadowing fading LoS
    X_sigma_NLoS = 0;       %% Shadowing fading NLoS

    %% Pathloss for LoS transmission
    PL_LoS = 32.4 + 21.6*log10(d)+20*log10(fc)+ X_sigma_LoS;

    %% Pathloss for NLoS transmission 
    PL_NLoS = 32.4 + 27.5*log10(d)+20*log10(fc)+ X_sigma_NLoS;
    
    %% Patloss function
    PL_CI = probability_LoS(d)*PL_LoS + (1-probability_LoS(d))*PL_NLoS;