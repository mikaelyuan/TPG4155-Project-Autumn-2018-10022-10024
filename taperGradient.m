function G = taperGradient(G)
    %% Taper gradient near sources to avoid interference
    [~,n] = size(G);
    taper = sin(linspace(0,pi/2,16));
    taper = repmat(taper,n,1)';   
    G(100:115,:) = taper.*G(100:115,:);
    
    %% Normalise gradient
    scale = 1.0/max(abs(G(:)));
    G = scale*G; 
end

