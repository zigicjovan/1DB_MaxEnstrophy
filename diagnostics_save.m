function diagnostics_save(J,K,step,momentum,timept,lambda,E0,ITER,testcase)
%{
VERBOSE 1: Add master diagnostics file and optimal initial conditions in physical space:
(1) [ Iterations, Enstrophy, Energy, Search Step Size, Search Momentum ]
<=>
[ n(iter), E_T(phi^(n)), K_T(phi^(n)), tau_n, beta_n ]
%}
    
    % Make master diagnostics file
    %[ n(iter), E_T(phi^(n)), K_T(phi^(n)), tau_n, beta_n ]
    diagnostics = zeros(ITER, 5); % make 5 column matrix
    diagnostics(:,1) = 1:ITER; % iteration number in col 1
    diagnostics(:,2) = J(1:ITER); % objective functional value in col 2
    diagnostics(:,3) = K(1:ITER); % energy in col 3
    diagnostics(:,4) = step(1:ITER); % search step in col 4
    diagnostics(:,5) = momentum(1:ITER); % search momentum in col 5

    master_file = [pwd '/data/diagnostics/diagnostics_E0_' num2str(E0)...
        '/diag' testcase '_E0(' num2str(E0) ')_lambda(' num2str(lambda) ').dat'];

    % append new columns
    try
        current = readmatrix(master_file);

        current_size = [ size(current,1) , size(diagnostics,1) ];
        size_choice = max(current_size);
        diag_new = NaN(size_choice, 5*timept); % make 5 column matrix
        
        if size(current,2) >= size(diag_new,2)
            diag_new = readmatrix(master_file);
        elseif size(current,2) < size(diag_new,2)
            diag_new( 1:size(current,1) , 1:((5*timept)-5) ) = readmatrix(master_file); 
        end

        % append new values
        diag_new( 1:size(diagnostics,1) , (5*timept)-4 ) = diagnostics(:,1); 
        diag_new( 1:size(diagnostics,1) , (5*timept)-3 ) = diagnostics(:,2); 
        diag_new( 1:size(diagnostics,1) , (5*timept)-2 ) = diagnostics(:,3); 
        diag_new( 1:size(diagnostics,1) , (5*timept)-1 ) = diagnostics(:,4); 
        diag_new( 1:size(diagnostics,1) , 5*timept ) = diagnostics(:,5); 
        diagnostics = diag_new;
    catch
        %
    end

    % Save files
    writematrix(diagnostics, master_file,'Delimiter','tab');

return