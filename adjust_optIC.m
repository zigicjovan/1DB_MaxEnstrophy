function result = adjust_optIC(phi,E,K0,N) 

    % Retraction to constraint manifold - see Ayala MS Thesis (3.23)
    
    % correction of non-zero mean solution
    if abs(mean(phi)) > 1e-5 
        phi(1) = 0;
    end

    phi_x = (2*pi*1i*K0).*phi;
    E0 = 0.5*sum(abs(phi_x).^2)/N^2;
    result = sqrt(E/E0)*phi;

return