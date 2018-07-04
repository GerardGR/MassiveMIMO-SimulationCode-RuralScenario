%% -------------------------------------------------------------
%% ------ Pathloss Representation: SNR as a function of d ------
%% -------------------------------------------------------------

tic;
%%
%% Pathloss Models: Sub-6GHz Frequencies
%%

distance_v = 10:10:5000;

figure('NumberTitle','off','Name','Pathloss Models: Sub-6GHz Frequencies');
 for ii=1:length(distance_v)
     distance = distance_v(ii);
     
     [Pathloss_vs_d_.Sub6GHz_.PL_M1(ii)] = pathloss_450MHz_6GHz_v1(distance, Conf_.MBS_ant_height, Conf_.UE_ant_height, ...
         Conf_.h, Conf_.W, Conf_.Freq_Carr);
     [Pathloss_vs_d_.Sub6GHz_.PL_CI(ii)] = pathloss_CI_v1(distance, Conf_.Freq_Carr);
     [Pathloss_vs_d_.Sub6GHz_.PL_CIH(ii)] = pathloss_CIH_v1(distance, Conf_.Freq_Carr, Conf_.MBS_ant_height);
 end
plot(distance_v,Pathloss_vs_d_.Sub6GHz_.PL_M1);
hold on
plot(distance_v,Pathloss_vs_d_.Sub6GHz_.PL_CI, '--'); 
plot(distance_v, Pathloss_vs_d_.Sub6GHz_.PL_CIH, '-.');
hold off

xlabel('d (m)');
ylabel('PL (dB)');
title('PathLoss Models - 5 GHz');
grid

%%
%% Pathloss Models: mmWave Frequency Bands
%%

figure('NumberTitle','off','Name','Pathloss Models: mmWave Frequency Bands');
 for ii=1:length(distance_v)
     distance = distance_v(ii);
     
     [Pathloss_vs_d_.mmW_.PL_CI(ii)] = pathloss_CI_v1(distance, Conf_.Freq_Carr2);
     [Pathloss_vs_d_.mmW_.PL_CIH(ii)] = pathloss_CIH_v1(distance, Conf_.Freq_Carr2, Conf_.MBS_ant_height);
 end

plot(distance_v,Pathloss_vs_d_.mmW_.PL_CI, '--'); 
hold on
plot(distance_v,Pathloss_vs_d_.mmW_.PL_CIH,'-.'); 
hold off

xlabel('d (m)');
ylabel('PL (dB)');
title('PathLoss Models - 28 GHz');
grid

toc;