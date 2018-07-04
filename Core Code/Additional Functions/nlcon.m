function [c,ceq] = nlcon(x,SNRmin, const_SNR)
    % alpha1: x(1), alpha2: x(2), gamma1: x(3), gamma2: x(4)

    SNR1 = (x(3)/x(1))*const_SNR;       %% Lineal
    SNR2 = (x(4)/x(2))*const_SNR;       %% Lineal
    
    c(1) = x(1)+ x(2)-1;
    c(2) = SNRmin - SNR1;
    c(3) = SNRmin - SNR2;
    c(4) = x(3)+ x(4)-1;
    ceq = [];
