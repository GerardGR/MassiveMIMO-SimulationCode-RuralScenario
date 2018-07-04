function [g] = g_function(sigma, phi, dH, M, K)
    if K > 1
        if sin(sigma) == sin(phi)
            g = M;
        else
            g = (sin(pi*dH*M*(sin(sigma)-sin(phi))).^2)/(M*sin(pi*dH*(sin(sigma)-sin(phi))).^2);
        end
    else 
        g = 0;
    end

