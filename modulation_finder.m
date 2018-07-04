function [M] = modulation_finder(x)
%     if (x/1024) > 1
%         M = 1024;
    if (x/256) > 1
        M = 256;
    elseif (x/64) > 1
        M = 64;
    elseif (x/16) > 1
        M = 16;
    elseif (x/4) > 1
        M = 4;
    else
        M = 0;
    end
        