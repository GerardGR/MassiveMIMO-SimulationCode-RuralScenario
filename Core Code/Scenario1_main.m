%% Scenario 1: Two users located at the same distance and with the same step
%% -- Rural Macro (RMa)
%% -- Analysis for the following carrier frequencies: 450 MHz, 5 GHz and 28 GHz
%% -- 

clear all;

%% ------------------------------------------------------------------------
%% Definition of the Scenario: Initial Assumptions
%% ------------------------------------------------------------------------
%% Pathloss and System Parameters
Conf_.MBS_ant_height=35;         %% m  antenna height
Conf_.UE_ant_height= 1.5;        %% m  antenna height

Conf_.NumUEs = 2;                        %% Number of users in the cell
Conf_.MBS_num_ant_tx = [1 25 50 100];    %% Number of Antennas at BS
Conf_.Mod = [4 16 64 256];               %% Vector of Modulations

Conf_.alpha_v = 0.01:0.01:0.99;     %% fraction of BW used
Conf_.gamma_v = 0.0:0.001:0.999;    %% fraction of Power used

Conf_.h = 5;
Conf_.W = 20;
Conf_.c = 3e8;
Conf_.d = 200;                   %% Distance of study for the Capacity Regions

Conf_.MBS_Pmax_dB =  46;         %% Max Power Transmitted by the BS
Conf_.UE_Pmax_dB =  23;          %% Max Power Transmitted by the UE
Conf_.BER = 1e-6;                %% Bit Error Rate

Conf_.BW = 5;                    %% Bandwidth, MHz [5 MHz, 20 MHz]
Conf_.BW2 = 100;                 %% Bandwidth, MHz [100 MHz]
Conf_.Freq_Carr= 0.45;           %% Carrier frequency, GHz [450 MHz, 5 GHz]
Conf_.Freq_Carr2= 28;            %% Carrier frequency, GHz [28 GHz]
Conf_.No = -174;                 %% Thermal noise den, dBm/Hz 
Conf_.F_UTnoise = 7;             %% Noise Figure, dB
Conf_.SNRmin = 10*log10(((4-1)/-1.6)*log(Conf_.BER/0.2)); %% SNRmin as a QPSK, dB

%% Massive MIMO: equation's variables
Conf_.dH = 0.5;     %% antenna spacing (in waveslengths)
Conf_.beta = 0;     %% describes the macroscopic large-scale fading. 0 <= beta <= 1
                    %%-> beta = 0, correspond to a negligibly weak inter-cell interference
                    %%-> beta = 1, means that the interference(inter-cell) is as strong as the desired signals.
%%
%% Noise
%%
Conf_.Eff_Noise  = Conf_.No + 10*log10(Conf_.BW*1e6)+Conf_.F_UTnoise;
Conf_.Eff_Noise2  = Conf_.No + 10*log10(Conf_.BW2*1e6)+Conf_.F_UTnoise;

%%
%% Path Loss Computation: 
%% 
    %%Path Loss: Model 1(450 MHz to 6 GHz)
    [Pathloss_.M1] = pathloss_450MHz_6GHz_v1(Conf_.d, Conf_.MBS_ant_height, Conf_.UE_ant_height, ...
                                             Conf_.h, Conf_.W, Conf_.Freq_Carr);
    %%Path Loss: CI (> 6 GHz)
    [Pathloss_.CI] = pathloss_CI_v1(Conf_.d, Conf_.Freq_Carr);
    %%Path Loss: CIH (> 6 GHz)
    [Pathloss_.CIH] = pathloss_CIH_v1(Conf_.d, Conf_.Freq_Carr, Conf_.MBS_ant_height);
    
%% ----------------------------------------
%% ------ CASE: Multiple Access TDMA ------
%% ----------------------------------------
%%
run Case_TDMA.m;    %%1-> Capacity Region
                    %%2-> SumRate as a function of d
                    %%3-> Modulation

%% ----------------------------------------
%% ------ CASE: Multiple Access SDMA ------
%% ----------------------------------------
%%
run Case_SDMA.m;    %%1-> Capacity Region
                    %%2-> Capacity Region 2 (SNRmin constraint)
                    %%3-> SumRate as a function of d
                    %%4-> Modulation
                    
%% -----------------------------------------------------------------
%% ------ CASE: Massive MIMO SDMA: (Intracell Interference) --------
%% -----------------------------------------------------------------
run MMIMO_SDMA.m;   %%1-> Evolution of the 'g' function given several sigmas
                    %%2-> MMIMO SDMA: Sum Rate LoS case as a function of d
                    
%% ----------------------------------------
%% ------ CASE: Multiple Access FDMA ------
%% ----------------------------------------
%%
run Case_FDMA.m;    %%1-> Capacity Region
                    %%2-> SumRate as a function of d
                    %%3-> Modulation
                    
%% -------------------------------------------------------------
%% ------ Pathloss Representation: SNR as a function of d ------
%% -------------------------------------------------------------
%%
run PathLoss_vs_d.m;    %%1-> Pathloss Models: Sub-6GHz Frequencies
                        %%2-> Pathloss Models: mmWave Frequency Bands

%% -------------------------------------------------
%% ------ Modulation Rate: as a function of d ------
%% -------------------------------------------------
%%
run Mod_evaluation_MA.m;  %%1-> Optimal Modulation evaluation for FDMA, TDMA and SDMA cases

%% -------------------------------
%% ------ Energy Efficiency ------
%% -------------------------------
%%
run EnergyEff.m;  %%1-> Effective Transmit Power
                  %%2-> Circuit Power
                  %%3-> Power Consumption
                  %%4-> Energy Efficiency: TDMA, SDMA and FDMA

