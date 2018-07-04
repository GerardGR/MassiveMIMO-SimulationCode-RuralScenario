%% ----------------------------------------
%% ------ CASE: Multiple Access TDMA ------
%% ----------------------------------------

tic;
%%
%% TDMA - Capacity Region
%%

TDMA_.Capacity_.R1 = zeros(1,length(Conf_.alpha_v));
TDMA_.Capacity_.R2 = zeros(1,length(Conf_.alpha_v));

TDMA_.Capacity_.SNR1 = Conf_.MBS_Pmax_dB-Pathloss_.CIH-Conf_.Eff_Noise;     % dB
TDMA_.Capacity_.SNR2 = Conf_.MBS_Pmax_dB-Pathloss_.CIH-Conf_.Eff_Noise;     % dB

if ((TDMA_.Capacity_.SNR1 >= Conf_.SNRmin) && (TDMA_.Capacity_.SNR2 >= Conf_.SNRmin))
    TDMA_.Capacity_.R1 = Conf_.alpha_v*Conf_.BW*log2(1+10^(TDMA_.Capacity_.SNR1/10));
    TDMA_.Capacity_.R2 = (1-Conf_.alpha_v)*Conf_.BW*log2(1+10^(TDMA_.Capacity_.SNR2/10));    
end

figure('NumberTitle','off','Name','TDMA: Capacity Region');
plot(TDMA_.Capacity_.R1,TDMA_.Capacity_.R2);
xlabel('R1');
ylabel('R2');
title('TDMA');
grid

%%
%% TDMA - SumRate as a function of d 
%%
distance_v = 10:10:5000;        %% vector of distance, step of 10 m

TDMA_.SumRate_vs_d_.SNR = zeros(length(Conf_.MBS_num_ant_tx),length(distance_v));
TDMA_.SumRate_vs_d_.Rate = zeros(length(Conf_.MBS_num_ant_tx),length(distance_v));

figure('NumberTitle','off','Name','TDMA: Rate as a function of d');
    
for jj=1:length(distance_v)
    distance = distance_v(jj);
    
    [TDMA_.SumRate_vs_d_.PL] = pathloss_CIH_v1(distance, Conf_.Freq_Carr, Conf_.MBS_ant_height);
    
    for ii=1:length(Conf_.MBS_num_ant_tx)
        M = Conf_.MBS_num_ant_tx(ii);
        
        TDMA_.SumRate_vs_d_.SNR(ii,jj) = M*10^((Conf_.MBS_Pmax_dB-TDMA_.SumRate_vs_d_.PL-Conf_.Eff_Noise)/10);
        if (10*log10(TDMA_.SumRate_vs_d_.SNR(ii,jj)) >= Conf_.SNRmin)
            TDMA_.SumRate_vs_d_.Rate(ii,jj) = Conf_.BW*log2(1+TDMA_.SumRate_vs_d_.SNR(ii,jj));  % Sum-rate (alpha1 + alpha2 = 1)
        else
            TDMA_.SumRate_vs_d_.Rate(ii,jj) = 0;
        end
    end
 end

plot(distance_v,TDMA_.SumRate_vs_d_.Rate);
xlabel('d');
ylabel('Sum-Rate');
legend(num2str(Conf_.MBS_num_ant_tx'));
title('TDMA - distance');
grid

%%
%% TDMA - Modulation
%%

TDMA_.Mod_.R1 = zeros(1,length(Conf_.alpha_v));  
TDMA_.Mod_.R2 = zeros(1,length(Conf_.alpha_v));

TDMA_.Mod_.const_SNR = Conf_.MBS_Pmax_dB-Pathloss_.CIH-Conf_.Eff_Noise;

figure('NumberTitle','off','Name','TDMA: Modulation');

for jj=1:length(Conf_.Mod)
    Mod2 = Conf_.Mod(jj); 
    for kk=1:length(Conf_.Mod)
        Mod1 = Conf_.Mod(kk);
        
        TDMA_.Mod_.SNR_mod1 = ((Mod1-1)/-1.6)*log(Conf_.BER/0.2);
        TDMA_.Mod_.SNR_mod2 = ((Mod2-1)/-1.6)*log(Conf_.BER/0.2);
        
        TDMA_.Mod_.R1_aux = Conf_.alpha_v*Conf_.BW*log2(Mod1);    
        TDMA_.Mod_.R2_aux = (1-Conf_.alpha_v)*Conf_.BW*log2(Mod2); 
              
        if ((10*log10(TDMA_.Mod_.SNR_mod1) >= Conf_.SNRmin) && (10*log10(TDMA_.Mod_.SNR_mod2) >= Conf_.SNRmin))
            if (TDMA_.Mod_.const_SNR >= 10*log10(TDMA_.Mod_.SNR_mod1)) && (TDMA_.Mod_.const_SNR >= 10*log10(TDMA_.Mod_.SNR_mod2))
                TDMA_.Mod_.R1 = TDMA_.Mod_.R1_aux;
                TDMA_.Mod_.R2 = TDMA_.Mod_.R2_aux;
            end
        end
        
        plot(TDMA_.Mod_.R1,TDMA_.Mod_.R2);
        hold on
    end
end

hold off
xlabel('R1');
ylabel('R2');
title('TDMA - Modulation');
grid

toc;