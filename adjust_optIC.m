function result = adjust_optIC(phi,E,K0,N) 
    
    phi_x = (2*pi*1i*K0).*phi;
    E0 = 0.5*sum(abs(phi_x).^2)/N^2;
    result = sqrt(E/E0)*phi;

return