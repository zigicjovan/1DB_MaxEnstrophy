function phi_save(phi, n, E0)
    
    len = length(phi);
    initialguess_file = [pwd '/optIC/phi_E0_' num2str(E0) '_' num2str(len) '_2048cont0.dat'];
    if n > 1
        previous_phi = readmatrix(initialguess_file);
        phi_initialguess = [ previous_phi ; phi ];
    else
        phi_initialguess = phi;
    end
    writematrix(phi_initialguess, initialguess_file,'Delimiter','tab');

return