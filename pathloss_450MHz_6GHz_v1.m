function [PL_1] = pathloss_450MHz_6GHz_v1(d, h_BS, h_UT, h, W, fc)
%%
%% distance_  in m
%% h_BS in m
%% h_UT in m
%% h, W in m
%% fc in GHz
%% derived from REP ITU-R M.2135-1, Table A1-2, pp 32
%% shadowing fading std  LoS (a): 4 dB
%%                       LoS (b): 6 dB                       
%%                       NLoS: 8 dB
%%

    %% Pathloss for LoS transmission
    d_BP = 2*pi*h_BS*h_UT*fc*1e9/3e8;
    if (d < d_BP) 
        %-- 10m < d < d_BP
        PL_LoS = 20*log10(40*pi*d*fc/3)+min(0.03*h^1.72,10)*log10(d)-min(0.044*h^1.72,14.77)+0.002*log10(h)*d;
    else
        %-- d_BP < d < 10000m; PL1(d_BP)
        PL1_BP = 20*log10(40*pi*d_BP*fc/3)+min(0.03*h^1.72,10)*log10(d_BP)-min(0.044*h^1.72,14.77)+0.002*log10(h)*d_BP;
        PL_LoS = PL1_BP + 40*log10(d/d_BP);
    end

    %% Pathloss for NLoS transmission
    %-- 10m < d < 5000m
    PL_NLoS = 161.04-7.1*log10(W)+7.5*log10(h)-(24.37-3.7*(h/h_BS)^2)*log10(h_BS)+(43.42-3.1*log10(h_BS))*(log10(d)-3)+20*log10(fc)-(3.2*(log10(11.75*h_UT))^2-4.97);
    
    %% Patloss function
    PL_1 = probability_LoS(d)*PL_LoS + (1-probability_LoS(d))*PL_NLoS;