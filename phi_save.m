function phi_save(phi, timept, E0, time_start, testcase)
    
    initialguess_file = [pwd '/optIC/phi_E0_' num2str(E0) '_' testcase '.dat'];
    if timept ~= time_start
        previous_phi = readmatrix(initialguess_file);
        phi_initialguess = [ previous_phi ; phi ];
    else
        phi_initialguess = phi;
    end
    writematrix(phi_initialguess, initialguess_file,'Delimiter','tab');

return