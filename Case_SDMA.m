%% ----------------------------------------
%% ------ CASE: Multiple Access SDMA ------
%% ----------------------------------------

tic;
%%
%% SDMA - Capacity Region 1
%%

SDMA_.Capacity_.Rmin = Conf_.BW*log2(1+10^(Conf_.SNRmin/10));
    
SDMA_.Capacity_.SNR1 = Conf_.gamma_v*10^((Conf_.MBS_Pmax_dB-Pathloss_.CIH-Conf_.Eff_Noise)/10);      %% Lineal
SDMA_.Capacity_.SNR2 = (1-Conf_.gamma_v)*10^((Conf_.MBS_Pmax_dB-Pathloss_.CIH-Conf_.Eff_Noise)/10);  %% Lineal

SDMA_.Capacity_.R1 = Conf_.BW*log2(1+SDMA_.Capacity_.SNR1);
SDMA_.Capacity_.R2 = Conf_.BW*log2(1+SDMA_.Capacity_.SNR2);

figure('NumberTitle','off','Name','SDMA: Capacity Region');
plot(SDMA_.Capacity_.R1,SDMA_.Capacity_.R2);
hold on

Q1 = [0 SDMA_.Capacity_.Rmin]; Q2 = [1000 SDMA_.Capacity_.Rmin];
plot([Q1(1) Q2(1)],[Q1(2) Q2(2)],'k', 'LineWidth',2);
hold on
Q3 = [SDMA_.Capacity_.Rmin 0]; Q4 = [SDMA_.Capacity_.Rmin 1000];
plot([Q3(1) Q4(1)],[Q3(2) Q4(2)],'k', 'LineWidth',2);
hold on

xlabel('R1');
ylabel('R2');
title('SDMA');
hold off
grid

%%
%% SDMA - Capacity Region 2 (SNRmin constraint)
%%

SDMA_.Capacity_constr_.R1 = zeros(1,length(Conf_.gamma_v));
SDMA_.Capacity_constr_.R2 = zeros(1,length(Conf_.gamma_v));

SNRmin_lineal = 10^(Conf_.SNRmin/10);  %% Lineal

for jj=1:length(Conf_.gamma_v)
    if (SDMA_.Capacity_.SNR1(jj) >= SNRmin_lineal) && (SDMA_.Capacity_.SNR2(jj) >= SNRmin_lineal)
        SDMA_.Capacity_constr_.R1(jj) = SDMA_.Capacity_.R1(jj);
        SDMA_.Capacity_constr_.R2(jj) = SDMA_.Capacity_.R2(jj);
    end
end
SDMA_.Capacity_constr_.R1(find(SDMA_.Capacity_constr_.R1==0)) = [];
SDMA_.Capacity_constr_.R2(find(SDMA_.Capacity_constr_.R2==0)) = [];

figure('NumberTitle','off','Name','SDMA: Capacity Region 2');
plot(SDMA_.Capacity_constr_.R1,SDMA_.Capacity_constr_.R2, 'r');
axis([0 1000 0 1000]);
xlabel('R1');
ylabel('R2');
title('SDMA');
grid

%%
%% SDMA - SumRate as a function of d 
%%
distance_v = 10:10:5000;        %% vector of distance, step of 10 m

%%Both users at the same distance -> optimal gamma = 0.5
SDMA_.SumRate_vs_d_.gamma_opt = 0.5;
SDMA_.SumRate_vs_d_.SNR1 = zeros(length(Conf_.MBS_num_ant_tx),length(distance_v));
SDMA_.SumRate_vs_d_.SNR2 = zeros(length(Conf_.MBS_num_ant_tx),length(distance_v));
SDMA_.SumRate_vs_d_.Rate = zeros(length(Conf_.MBS_num_ant_tx),length(distance_v));

figure('NumberTitle','off','Name','SDMA: Rate as a function of d');

for jj=1:length(distance_v)
    distance = distance_v(jj);
    
    [SDMA_.SumRate_vs_d_.PL] = pathloss_CIH_v1(distance, Conf_.Freq_Carr, Conf_.MBS_ant_height);
        
    for ii=1:length(Conf_.MBS_num_ant_tx)
        M = Conf_.MBS_num_ant_tx(ii);
        
        SDMA_.SumRate_vs_d_.SNR1(ii,jj) = M*SDMA_.SumRate_vs_d_.gamma_opt...
                                          *10^((Conf_.MBS_Pmax_dB-SDMA_.SumRate_vs_d_.PL-Conf_.Eff_Noise)/10); 
        SDMA_.SumRate_vs_d_.SNR2(ii,jj) = M*(1-SDMA_.SumRate_vs_d_.gamma_opt)...
                                          *10^((Conf_.MBS_Pmax_dB-SDMA_.SumRate_vs_d_.PL-Conf_.Eff_Noise)/10);     

        if (SDMA_.SumRate_vs_d_.SNR1(ii,jj) >= SNRmin_lineal) && (SDMA_.SumRate_vs_d_.SNR2(ii,jj) >= SNRmin_lineal)
            SDMA_.SumRate_vs_d_.Rate(ii,jj) = Conf_.BW*log2(1+SDMA_.SumRate_vs_d_.SNR1(ii,jj)) + ...
                                              Conf_.BW*log2(1+SDMA_.SumRate_vs_d_.SNR2(ii,jj));  
        else
            SDMA_.SumRate_vs_d_.Rate(ii,jj) = 0;
        end
    end
end

plot(distance_v,SDMA_.SumRate_vs_d_.Rate);
xlabel('d');
ylabel('Sum-Rate');
legend(num2str(Conf_.MBS_num_ant_tx'));
title('SDMA - distance');
grid

%%
%% SDMA - Modulation
%%

SDMA_.Mod_.R1 = zeros(1,length(Conf_.gamma_v));
SDMA_.Mod_.R2 = zeros(1,length(Conf_.gamma_v));

figure('NumberTitle','off','Name','SDMA: Modulation');

for jj=1:length(Conf_.Mod)
    Mod2 = Conf_.Mod(jj); 
    for kk=1:length(Conf_.Mod)
        Mod1 = Conf_.Mod(kk);
        SDMA_.Mod_.SNR_mod1 = ((Mod1-1)/-1.6)*log(Conf_.BER/0.2);
        SDMA_.Mod_.SNR_mod2 = ((Mod2-1)/-1.6)*log(Conf_.BER/0.2);
        
        SDMA_.Mod_.R1_aux = Conf_.BW*log2(Mod1); %modulació M1 and uncoded
        SDMA_.Mod_.R2_aux = Conf_.BW*log2(Mod2);
        
        SDMA_.Mod_.P1_aux = 10^((10*log10(SDMA_.Mod_.SNR_mod1)+Pathloss_.CIH+Conf_.Eff_Noise)/10);
        SDMA_.Mod_.P2_aux = 10^((10*log10(SDMA_.Mod_.SNR_mod2)+Pathloss_.CIH+Conf_.Eff_Noise)/10);
        
        if (SDMA_.Mod_.P1_aux + SDMA_.Mod_.P2_aux <= 10^(Conf_.MBS_Pmax_dB/10)) 
            SDMA_.Mod_.R1 = SDMA_.Mod_.R1_aux;
            SDMA_.Mod_.R2 = SDMA_.Mod_.R2_aux;
        end
        scatter(SDMA_.Mod_.R1,SDMA_.Mod_.R2);
        hold on
    end
end
hold off
xlabel('R1');
ylabel('R2');
title('SDMA - Modulation');
axis([0 150 0 150]);
grid

toc;