%% -----------------------------------------------------------
%% ------ Modulation Rate: as a function of d ----------------
%% -----------------------------------------------------------

tic;
%%
%% Optimal Modulation evaluation for FDMA, TDMA and SDMA cases
%%

distance_v = 10:10:5000;    %% distance of interest

 for jj=1:length(distance_v)
     distance = distance_v(jj);
     
     [Mod_test_.PL_M1(jj)] = pathloss_450MHz_6GHz_v1(distance, Conf_.MBS_ant_height, Conf_.UE_ant_height, ...
         Conf_.h, Conf_.W, Conf_.Freq_Carr);
     [Mod_test_.PL_CI(jj)] = pathloss_CI_v1(distance, Conf_.Freq_Carr);
     [Mod_test_.PL_CIH(jj)] = pathloss_CIH_v1(distance, Conf_.Freq_Carr, Conf_.MBS_ant_height);
     
     for ii=1:length(Conf_.MBS_num_ant_tx)
         M = Conf_.MBS_num_ant_tx(ii);
        
         Mod_test_.FDMA_.SNR1(ii,jj) = (FDMA_.SumRate_vs_d_.gamma1_opt(ii,jj)/FDMA_.SumRate_vs_d_.alfa1_opt(ii,jj))*...
            M*10^(0.1*(Conf_.MBS_Pmax_dB-Mod_test_.PL_CIH(jj)-Conf_.Eff_Noise));

         Mod_test_.FDMA_.SNR2(ii,jj) = (FDMA_.SumRate_vs_d_.gamma2_opt(ii,jj)/FDMA_.SumRate_vs_d_.alfa2_opt(ii,jj))*...
             M*10^(0.1*(Conf_.MBS_Pmax_dB-Mod_test_.PL_CIH(jj)-Conf_.Eff_Noise));
     
         Mod_test_.SDMA_.SNR(ii,jj) = M*SDMA_.SumRate_vs_d_.gamma_opt* ...
             10^((Conf_.MBS_Pmax_dB-Mod_test_.PL_CIH(jj)-Conf_.Eff_Noise)/10);
         
         Mod_test_.TDMA_.SNR(ii,jj) = M*10^((Conf_.MBS_Pmax_dB-Mod_test_.PL_CIH(jj)-Conf_.Eff_Noise)/10);
     
         Mod_test_.FDMA_.Mod_aux1(ii,jj) = (-1.6*Mod_test_.FDMA_.SNR1(ii,jj))/(log(Conf_.BER/0.2))+1;
         Mod_test_.FDMA_.Mod_aux2(ii,jj) = (-1.6*Mod_test_.FDMA_.SNR2(ii,jj))/(log(Conf_.BER/0.2))+1;
         Mod_test_.SDMA_.Mod_aux(ii,jj) = (-1.6*Mod_test_.SDMA_.SNR(ii,jj))/(log(Conf_.BER/0.2))+1;
         Mod_test_.TDMA_.Mod_aux(ii,jj) = (-1.6*Mod_test_.TDMA_.SNR(ii,jj))/(log(Conf_.BER/0.2))+1;
         
         Mod_test_.FDMA_.Modulation1(ii,jj) = modulation_finder(Mod_test_.FDMA_.Mod_aux1(ii,jj));
         Mod_test_.FDMA_.Modulation2(ii,jj) = modulation_finder(Mod_test_.FDMA_.Mod_aux2(ii,jj));
         Mod_test_.SDMA_.Modulation(ii,jj) = modulation_finder(Mod_test_.SDMA_.Mod_aux(ii,jj));
         Mod_test_.TDMA_.Modulation(ii,jj) = modulation_finder(Mod_test_.TDMA_.Mod_aux(ii,jj));
         
         Mod_test_.FDMA_.SumRate(ii,jj) = FDMA_.SumRate_vs_d_.alfa1_opt(ii,jj)*Conf_.BW* ...
             log2(Mod_test_.FDMA_.Modulation1(ii,jj))+ FDMA_.SumRate_vs_d_.alfa2_opt(ii,jj)* ...
             Conf_.BW*log2(Mod_test_.FDMA_.Modulation2(ii,jj));
         Mod_test_.SDMA_.SumRate(ii,jj) = 2*Conf_.BW*log2(Mod_test_.SDMA_.Modulation(ii,jj));
         Mod_test_.TDMA_.SumRate(ii,jj) = Conf_.BW*log2(Mod_test_.TDMA_.Modulation(ii,jj));
     end
 end

toc;