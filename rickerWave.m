function [ricker] = rickerWave(freq,dims)
%RICKERWAVE Generates the time derivative of a Ricker wavelet
% to be used as a source signature
    t = 0:dims.dt:dims.dt*(dims.nt);
    tau  = (t-1/freq)*freq*pi;

    ricker = (1-tau.*tau).*exp(-tau.*tau);
    ricker = (ricker(2:end) - ricker(1:end-1))';
end

