function result = adjust_optIC(phi,E,K0,N) 
    
    if abs(mean(phi)) > 1e-5 % correction of non-zero mean solution
        phi(1) = 0;
        %{
        phi_phys = real(ifft(phi));
        shift_ic = phi_phys(N/2);
        phi_phys = phi_phys - shift_ic;
        firsthalf = abs(sum(phi_phys(1:N/2)));
        secondhalf = abs(sum(phi_phys((N/2+1):end)));
        scalefactor = firsthalf/secondhalf;
        phi_phys(1:N/2) = phi_phys(1:N/2) / scalefactor ;
        % phi_phys((N/2)+1:N) = phi_phys((N/2)+1:N) * scalefactor ;
        phi = fft(phi_phys);
        %}
    end
    phi_x = (2*pi*1i*K0).*phi;
    E0 = 0.5*sum(abs(phi_x).^2)/N^2;
    result = sqrt(E/E0)*phi;

return