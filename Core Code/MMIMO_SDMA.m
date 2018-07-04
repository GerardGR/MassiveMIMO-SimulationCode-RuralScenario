%% -----------------------------------------------------------------
%% ------ CASE: Massive MIMO SDMA: (Intracell Interference) --------
%% -----------------------------------------------------------------

tic;
%%
%% Evolution of the 'g' function given several sigmas
%%

MMIMO_SDMA_.G_function_.phi_v = 0:0.01:2*pi;    %% Angle of interfering UE 
MMIMO_SDMA_.G_function_.sigma_ref = pi/6;       %% Reference angle: 30º
MMIMO_SDMA_.G_function_.sigma_v = zeros(1,21);  %% Angle of interfering UE

G_M1 = zeros(length(MMIMO_SDMA_.G_function_.sigma_v),length(MMIMO_SDMA_.G_function_.phi_v));
G_M2 = zeros(length(MMIMO_SDMA_.G_function_.sigma_v),length(MMIMO_SDMA_.G_function_.phi_v));
G_M3 = zeros(length(MMIMO_SDMA_.G_function_.sigma_v),length(MMIMO_SDMA_.G_function_.phi_v));
G_M4 = zeros(length(MMIMO_SDMA_.G_function_.sigma_v),length(MMIMO_SDMA_.G_function_.phi_v));

for kk = 1:length(MMIMO_SDMA_.G_function_.sigma_v)
    if (kk == 1)
        MMIMO_SDMA_.G_function_.sigma_v(kk) = MMIMO_SDMA_.G_function_.sigma_ref;
    else
        x = rand;
        MMIMO_SDMA_.G_function_.sigma_v(kk)= 2*pi*x;
    end
end

%% Macro G function definition: for different antennas and the sigma-phi angles
MMIMO_SDMA_.G_function_.G = cat(3,G_M1,G_M2,G_M3,G_M4);

for kk=1:length(Conf_.MBS_num_ant_tx)
    M = Conf_.MBS_num_ant_tx(kk);
    for jj= 1:length(MMIMO_SDMA_.G_function_.sigma_v)
        sigma = MMIMO_SDMA_.G_function_.sigma_v(jj);
        for ii=1:length(MMIMO_SDMA_.G_function_.phi_v)
            phi = MMIMO_SDMA_.G_function_.phi_v(ii);
            MMIMO_SDMA_.G_function_.G(jj,ii,kk) = g_function(sigma, phi, Conf_.dH, M, Conf_.NumUEs);
        end
    end
end

   
%%
%% SDMA: Sum Rate LoS case as a function of d
%%
%% Achivable DL sum SE [bit/s/Hz/cell] in the LoS case

distance_v = 10:10:5000;        %% vector of distance, step of 10 m
SNRmin_lineal = 10^(Conf_.SNRmin/10);  %% Lineal

MMIMO_SDMA_.SumRate_vs_d_.gamma_opt = 0.5;
MMIMO_SDMA_.SumRate_vs_d_.eSum_SE_LoS = zeros(length(MMIMO_SDMA_.G_function_.sigma_v),length(MMIMO_SDMA_.G_function_.phi_v));
MMIMO_SDMA_.SumRate_vs_d_.Avg_eSum_SE_LoS = zeros(length(Conf_.MBS_num_ant_tx),length(distance_v));
MMIMO_SDMA_.SumRate_vs_d_.Avgg_eSum_SE_LoS = zeros(length(Conf_.MBS_num_ant_tx),length(distance_v));
MMIMO_SDMA_.SumRate_vs_d_.Avgl_eSum_SE_LoS = zeros(length(Conf_.MBS_num_ant_tx),length(distance_v));

Mod_rate_SE = zeros(length(MMIMO_SDMA_.G_function_.sigma_v),length(MMIMO_SDMA_.G_function_.phi_v));

bar = waitbar(0,'- Computing Massive MIMO SDMA: Avg Sum Rate LoS case -');

for jj=1:length(distance_v)
    waitbar(jj/length(distance_v),bar,'- Computing Massive MIMO SDMA: Avg Sum Rate LoS case -');
    distance = distance_v(jj);
    
    [MMIMO_SDMA_.SumRate_vs_d_.PL] = pathloss_CIH_v1(distance, Conf_.Freq_Carr, Conf_.MBS_ant_height);
    
    for ii=1:length(Conf_.MBS_num_ant_tx)
        M = Conf_.MBS_num_ant_tx(ii);
        
        MMIMO_SDMA_.SumRate_vs_d_.const_SNR = 10^((Conf_.MBS_Pmax_dB-MMIMO_SDMA_.SumRate_vs_d_.PL ...
            -Conf_.Eff_Noise)/10)/Conf_.NumUEs; 
        
        for kk= 1:length(MMIMO_SDMA_.G_function_.sigma_v)
            sigma = MMIMO_SDMA_.G_function_.sigma_v(kk);
            for nn=1:length(MMIMO_SDMA_.G_function_.phi_v)
                phi = MMIMO_SDMA_.G_function_.phi_v(nn);
                
                MMIMO_SDMA_.SumRate_vs_d_.SNR = MMIMO_SDMA_.SumRate_vs_d_.gamma_opt*M/...
                    (MMIMO_SDMA_.G_function_.G(kk,nn,ii)+ 1/MMIMO_SDMA_.SumRate_vs_d_.const_SNR);
                
                if (MMIMO_SDMA_.SumRate_vs_d_.SNR >= SNRmin_lineal)
                    MMIMO_SDMA_.SumRate_vs_d_.eSum_SE_LoS(kk,nn) = Conf_.BW*Conf_.NumUEs*log2(1+MMIMO_SDMA_.SumRate_vs_d_.SNR);
                    
                    Mod_aux = (-1.6*MMIMO_SDMA_.SumRate_vs_d_.SNR)/(log(Conf_.BER/0.2))+1;
                    Modulation = modulation_finder(Mod_aux);
                    Mod_rate_SE(kk,nn) = Conf_.NumUEs*Conf_.BW*log2(Modulation);
                else
                    MMIMO_SDMA_.SumRate_vs_d_.eSum_SE_LoS(kk,nn) =  0;
                    
                    Mod_rate_SE(kk,nn) = 0;
                end
            end
        end      
        mean_SE_max = mean(MMIMO_SDMA_.SumRate_vs_d_.eSum_SE_LoS(:));
        mean_SE_mod = mean(Mod_rate_SE(:));
        Mod_test_.MMIMO_SDMA_.SumRate(ii,jj) = mean_SE_mod;
                
        MMIMO_SDMA_.SumRate_vs_d_.Avg_eSum_SE_LoS(ii,jj) = mean_SE_max;
        
        i_greater = find(MMIMO_SDMA_.SumRate_vs_d_.eSum_SE_LoS - mean_SE_max > 0);
        i_lower = find(MMIMO_SDMA_.SumRate_vs_d_.eSum_SE_LoS - mean_SE_max < 0);
        
        m_greater = mean(MMIMO_SDMA_.SumRate_vs_d_.eSum_SE_LoS(i_greater));
        m_lower = mean(MMIMO_SDMA_.SumRate_vs_d_.eSum_SE_LoS(i_lower));
        
        var_greater = var(MMIMO_SDMA_.SumRate_vs_d_.eSum_SE_LoS(i_greater));
        var_lower = var(MMIMO_SDMA_.SumRate_vs_d_.eSum_SE_LoS(i_lower));
        
        MMIMO_SDMA_.SumRate_vs_d_.Avgg_eSum_SE_LoS(ii,jj) = sqrt(var_greater + m_greater^2);
        MMIMO_SDMA_.SumRate_vs_d_.Avgl_eSum_SE_LoS(ii,jj) = sqrt(var_lower + m_lower^2);
        
    end
end
close(bar);

figure('NumberTitle','off','Name','Massive MIMO SDMA: Avg Sum Rate LoS case');
plot(distance_v,MMIMO_SDMA_.SumRate_vs_d_.Avg_eSum_SE_LoS(2,:));
hold on
plot(distance_v,MMIMO_SDMA_.SumRate_vs_d_.Avgl_eSum_SE_LoS(2,:), ':');
plot(distance_v,MMIMO_SDMA_.SumRate_vs_d_.Avgg_eSum_SE_LoS(2,:), '--');
hold off
xlabel('d');
ylabel('Rate');
title('SDMA: Avg Sum Rate LoS case - distance');
grid

toc;
