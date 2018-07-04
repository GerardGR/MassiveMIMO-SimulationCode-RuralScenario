%% -------------------------------
%% ------ Energy Efficiency ------
%% -------------------------------
%%

tic;
%%
%% Effective Transmit Power
%%
EnergyEff_.mu_BS = 0.6;     %% efficiency of the power amplifier
EnergyEff_.mu_UE = 0.6;     %% efficiency of the power amplifier to the UE

EnergyEff_.ETP = 10^(0.1*(Conf_.MBS_Pmax_dB-30))/EnergyEff_.mu_BS + Conf_.NumUEs*10^(0.1*(Conf_.UE_Pmax_dB-30))/EnergyEff_.mu_UE;

%%
%% Circuit Power
%%
EnergyEff_.P_fix = 10;      %% constant value (10W)
EnergyEff_.PBS = 1;         %% Power to let the components run (1W)

EnergyEff_.CP = EnergyEff_.P_fix  + Conf_.MBS_num_ant_tx*EnergyEff_.PBS/10e6*Conf_.BW*1e6;

%%
%% Power Consumption
%%
EnergyEff_.PowerCons = EnergyEff_.ETP + EnergyEff_.CP;

%%
%% Energy Efficiency: TDMA, SDMA and FDMA
%%
%% Recalculating the CP and PC for FDMA due to BW allocation
EnergyEff_.FDMA_CP = EnergyEff_.P_fix  + Conf_.MBS_num_ant_tx'*EnergyEff_.PBS/10e6* ...
   Conf_.BW.*(FDMA_.SumRate_vs_d_.alfa1_opt + FDMA_.SumRate_vs_d_.alfa2_opt)*1e6;
EnergyEff_.FDMA_PowerCons = EnergyEff_.ETP + EnergyEff_.FDMA_CP;

%% EE at d using M antennas: maximum rate
EnergyEff_.TDMA_d = TDMA_.SumRate_vs_d_.Rate./EnergyEff_.PowerCons';
EnergyEff_.SDMA_d = SDMA_.SumRate_vs_d_.Rate./EnergyEff_.PowerCons';
EnergyEff_.MMIMO_SDMA_d = MMIMO_SDMA_.SumRate_vs_d_.Avg_eSum_SE_LoS./EnergyEff_.PowerCons';
EnergyEff_.FDMA_d = FDMA_.SumRate_vs_d_.Rate./EnergyEff_.FDMA_PowerCons;

%% EE at d using M antennas: modulation rate
EnergyEff_.TDMA_d_mod = Mod_test_.TDMA_.SumRate./EnergyEff_.PowerCons';
EnergyEff_.SDMA_d_mod = Mod_test_.SDMA_.SumRate./EnergyEff_.PowerCons';
EnergyEff_.MMIMO_SDMA_d_mod = Mod_test_.MMIMO_SDMA_.SumRate./EnergyEff_.PowerCons';
EnergyEff_.FDMA_d_mod = Mod_test_.FDMA_.SumRate./EnergyEff_.FDMA_PowerCons;

toc;