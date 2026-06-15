function phi_save(phi, timept, E0, testcase)                %
    
    initialguess_file = [pwd '/data/optIC/phi' testcase '_E0(' num2str(E0) ').dat'];     %
    if timept ~= 1                                 %
        previous_phi = readmatrix(initialguess_file);       %
        phi_initialguess = [ previous_phi ; phi ];          %
    else                            %
        phi_initialguess = phi;     %
    end                                                                     %
    writematrix(phi_initialguess, initialguess_file,'Delimiter','tab');     %

return