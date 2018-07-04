%% ----------------------------------------
%% ------ CASE: Multiple Access FDMA ------
%% ----------------------------------------

tic;
warning off;

%%
%% FDMA - Capacity Region
%%
%% Fmincon: Parameters' Initialization

Fmincon_.mu_1_v = 0:0.005:0.999;       %% optimitzation parameter
Fmincon_.lb = 0.01*ones(4);
Fmincon_.ub = 0.99*ones(4);
Fmincon_.A = [];
Fmincon_.b = [];
Fmincon_.Aeq = [];
Fmincon_.beq = [];
Fmincon_.x0 = [0.02, 0.02, 0.02, 0.02];
Fmincon_.options=optimoptions('fmincon');
Fmincon_.options.Display='off';

FDMA_.Capacity_.const_SNR = 10^((Conf_.MBS_Pmax_dB-Pathloss_.CIH-Conf_.Eff_Noise)/10);   %% Lineal
SNRmin_lineal = 10^(Conf_.SNRmin/10);   %% Lineal

FDMA_.Capacity_.R1 = zeros(1,length(Fmincon_.mu_1_v));
FDMA_.Capacity_.R2 = zeros(1,length(Fmincon_.mu_1_v));

figure('NumberTitle','off','Name','FDMA: Capacity Region');

for jj=1:length(Fmincon_.mu_1_v)
    mu_1=Fmincon_.mu_1_v(jj);
    
    % alpha1: x(1), alpha2: x(2), gamma1: x(3), gamma2: x(4)
    objective = @(x)(-mu_1*(x(1)*Conf_.BW*log2(1+(x(3)/x(1))*FDMA_.Capacity_.const_SNR)) ...
                -(1-mu_1)*x(2)*Conf_.BW*log2(1+(x(4)/x(2))*FDMA_.Capacity_.const_SNR));
            
    [x,fval,exitflag] = fmincon(objective,Fmincon_.x0,Fmincon_.A,Fmincon_.b,Fmincon_.Aeq, Fmincon_.beq, ...
        Fmincon_.lb,Fmincon_.ub, @(x)nlcon(x, SNRmin_lineal, FDMA_.Capacity_.const_SNR), Fmincon_.options);

    FDMA_.Capacity_.R1(jj)= x(1)*Conf_.BW*log2(1+(x(3)/x(1))*FDMA_.Capacity_.const_SNR);
    FDMA_.Capacity_.R2(jj)= x(2)*Conf_.BW*log2(1+(x(4)/x(2))*FDMA_.Capacity_.const_SNR);
end

plot(FDMA_.Capacity_.R1,FDMA_.Capacity_.R2,'x-');
xlabel('R1');
ylabel('R2');
title('FDMA - Max Sum Rate');
grid

%%
%% FDMA - SumRate as a function of d 
%%
%% Looking for the optimal alfa/gamma's vector

distance_v = 10:10:5000;        %% vector of distance, step of 10 m

FDMA_.SumRate_vs_d_.alfa1_opt = zeros(length(Conf_.MBS_num_ant_tx),length(distance_v));
FDMA_.SumRate_vs_d_.alfa2_opt = zeros(length(Conf_.MBS_num_ant_tx),length(distance_v));
FDMA_.SumRate_vs_d_.gamma1_opt = zeros(length(Conf_.MBS_num_ant_tx),length(distance_v));
FDMA_.SumRate_vs_d_.gamma2_opt = zeros(length(Conf_.MBS_num_ant_tx),length(distance_v));

FDMA_.SumRate_vs_d_.const_SNR = zeros(length(Conf_.MBS_num_ant_tx),length(distance_v));
FDMA_.SumRate_vs_d_.Rate = zeros(length(Conf_.MBS_num_ant_tx),length(distance_v));

R1_aux = zeros(1,length(Fmincon_.mu_1_v));
R2_aux = zeros(1,length(Fmincon_.mu_1_v));

bar_v2 = waitbar(0,'- FDMA: Rate as a function of d -');
for jj=1:length(distance_v)
    waitbar(jj/length(distance_v),bar_v2,'- FDMA: Rate as a function of d -');
    distance = distance_v(jj);
    
%     [FDMA_.SumRate_vs_d_.PL] = pathloss_450MHz_6GHz_v1(distance, Conf_.MBS_ant_height, ...
%         Conf_.UE_ant_height, Conf_.h, Conf_.W, Conf_.Freq_Carr);
    [FDMA_.SumRate_vs_d_.PL] = pathloss_CIH_v1(distance, Conf_.Freq_Carr, Conf_.MBS_ant_height);
    
    for ii=1:length(Conf_.MBS_num_ant_tx)
        M = Conf_.MBS_num_ant_tx(ii);
        FDMA_.SumRate_vs_d_.const_SNR(ii,jj) = M*10^((Conf_.MBS_Pmax_dB-FDMA_.SumRate_vs_d_.PL-Conf_.Eff_Noise)/10);
        Rmax = 0;
        
%% Organization of the optimization variables
%% First   NumUEs positions (1, NumUEs)  ------------> fraction of BW assigned to each UEs
%% Second  NumUEs positions (NumUEs+1, 2*NumUEs) ----> fraction of  Power allocated to each UE

        objective = @(x) -sum(x(1:Conf_.NumUEs)*Conf_.BW.*log2(1+(x(Conf_.NumUEs+1:2*Conf_.NumUEs)./ ...
            x(1:Conf_.NumUEs)).*FDMA_.SumRate_vs_d_.const_SNR(ii,jj).'));       

        [x,fval,exitflag] = fmincon(objective,Fmincon_.x0,Fmincon_.A,Fmincon_.b,Fmincon_.Aeq, Fmincon_.beq, ...
            Fmincon_.lb,Fmincon_.ub, @(x)nlcon(x, SNRmin_lineal, FDMA_.SumRate_vs_d_.const_SNR(ii,jj)),Fmincon_.options);

        FDMA_.SumRate_vs_d_.alfa1_opt(ii,jj) = x(1);    
        FDMA_.SumRate_vs_d_.alfa2_opt(ii,jj) = x(2);
        FDMA_.SumRate_vs_d_.gamma1_opt(ii,jj) = x(3);
        FDMA_.SumRate_vs_d_.gamma2_opt(ii,jj) = x(4);

%% Computing Sum-Rate using the optimal parameters
        FDMA_.SumRate_vs_d_.Rate(ii,jj) = FDMA_.SumRate_vs_d_.alfa1_opt(ii,jj)*Conf_.BW*log2(1+ ...
            (FDMA_.SumRate_vs_d_.gamma1_opt(ii,jj)/FDMA_.SumRate_vs_d_.alfa1_opt(ii,jj)) ...
            *FDMA_.SumRate_vs_d_.const_SNR(ii,jj)) + FDMA_.SumRate_vs_d_.alfa2_opt(ii,jj)*Conf_.BW...
            *log2(1+(FDMA_.SumRate_vs_d_.gamma2_opt(ii,jj)/FDMA_.SumRate_vs_d_.alfa2_opt(ii,jj))...
            *FDMA_.SumRate_vs_d_.const_SNR(ii,jj));

    end
end
close(bar_v2);

figure('NumberTitle','off','Name','FDMA: Rate as a function of d');
plot(distance_v,FDMA_.SumRate_vs_d_.Rate);
xlabel('d (m) ');
ylabel('Sum-Rate (Mbps)');
legend(num2str(Conf_.MBS_num_ant_tx'));
title('FDMA - distance');
grid

FDMA_.SumRate_vs_d_.BW1 = FDMA_.SumRate_vs_d_.alfa1_opt*Conf_.BW;
FDMA_.SumRate_vs_d_.BW2 = FDMA_.SumRate_vs_d_.alfa2_opt*Conf_.BW;

figure('NumberTitle','off','Name','FDMA: Optimal Parameters as a function of d');
ax1 = subplot(1,2,1); 
plot(ax1, distance_v,FDMA_.SumRate_vs_d_.BW1(1,:), distance_v, FDMA_.SumRate_vs_d_.BW2(1,:));
xlabel(ax1,'d');
ylabel(ax1,'Bandwidth');
title(ax1,'FDMA: Bandwidth - distance');
grid

ax2 = subplot(1,2,2); 
plot(ax2, distance_v,FDMA_.SumRate_vs_d_.gamma1_opt(1,:), distance_v, FDMA_.SumRate_vs_d_.gamma2_opt(1,:));
xlabel(ax2,'d');
ylabel(ax2,'Percentage of Power Allocation');
title(ax2,'FDMA: Power Allocation(%) - distance');
grid

%%
%% FDMA - Modulation
%%

FDMA_.Mod_.R1 = zeros(1,length(Fmincon_.mu_1_v));  
FDMA_.Mod_.R2 = zeros(1,length(Fmincon_.mu_1_v));

FDMA_.Mod_.const_SNR = 10^((Conf_.MBS_Pmax_dB-Pathloss_.CIH-Conf_.Eff_Noise)/10);

figure('NumberTitle','off','Name','FDMA: Modulation');

for jj=1:length(Conf_.Mod)
    Mod2 = Conf_.Mod(jj); 
    for ii=1:length(Conf_.Mod)
        Mod1 = Conf_.Mod(ii);
        
        for kk=1:length(Fmincon_.mu_1_v)
            mu_1=Fmincon_.mu_1_v(kk);
                     
            FDMA_.Mod_.SNR_mod1 = ((Mod1-1)/-1.6)*log(Conf_.BER/0.2);
            FDMA_.Mod_.SNR_mod2 = ((Mod2-1)/-1.6)*log(Conf_.BER/0.2);
            
            % alpha1: x(1), alpha2: x(2), gamma1: x(3), gamma2: x(4)
            objective = @(x)(-mu_1*(x(1)*Conf_.BW*log2(1+(x(3)/x(1))*FDMA_.Mod_.const_SNR)) ...
                        -(1-mu_1)*x(2)*Conf_.BW*log2(1+(x(4)/x(2))*FDMA_.Mod_.const_SNR));
            
            [x,fval,exitflag] = fmincon(objective,Fmincon_.x0,Fmincon_.A,Fmincon_.b,Fmincon_.Aeq, Fmincon_.beq, ...
                                Fmincon_.lb,Fmincon_.ub, @(x)nlcon2(x,FDMA_.Mod_.SNR_mod1, FDMA_.Mod_.SNR_mod2, ...
                                SNRmin_lineal, FDMA_.Mod_.const_SNR),Fmincon_.options);
                    
            FDMA_.Mod_.R1(kk) = x(1)*Conf_.BW*log2(Mod1);    
            FDMA_.Mod_.R2(kk) = x(2)*Conf_.BW*log2(Mod2); 

        end
        plot(FDMA_.Mod_.R1,FDMA_.Mod_.R2,'*-');
        hold on
    end
end
xlabel('R1');
ylabel('R2');
title('FDMA - Max Sum Rate - Modulation');
hold off
grid

toc;