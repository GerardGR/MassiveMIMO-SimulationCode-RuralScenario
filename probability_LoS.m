function [Prob_LoS] = probability_LoS(d)
    if (d <= 10)
        Prob_LoS = 1;
    else
        Prob_LoS = exp(-(d-10)/1000);
    end