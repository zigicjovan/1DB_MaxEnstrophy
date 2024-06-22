function diagnostics_save(J,K,step,momentum,timept,lambda,E0,ITER,testcase)
%{
VERBOSE 1: Add master diagnostics file and optimal initial conditions in physical space:
(1) [ Iterations, Enstrophy, Energy, Search Step Size, Search Momentum ]
<=>
[ n(iter), E_T(phi^(n)), K_T(phi^(n)), tau_n, beta_n ]
%}

    master_file = [pwd '/data/diagnostics/diag' testcase '_E0(' num2str(E0) ')_lambda(' num2str(lambda) ').dat'];

    % append new columns
    try
        current = readmatrix(master_file);
        current_size = [ size(current,1) , ITER ];
        size_choice = max(current_size);
        diag_new = NaN(size_choice, 5*timept); % make 5 column matrix
        diag_new( 1:size(current,1) , 1:size(current,2) ) = current; 
    catch
        diag_new = NaN(ITER, 5*timept); % make 2 column matrix
    end

    % append new values
    diag_new( 1:ITER , (5*timept)-4 ) = 1:ITER; % iteration number in col 1 
    diag_new( 1:ITER , (5*timept)-3 ) = J(1:ITER); % objective functional value in col 2 
    diag_new( 1:ITER , (5*timept)-2 ) = K(1:ITER); % energy in col 3
    diag_new( 1:ITER , (5*timept)-1 ) = step(1:ITER); % search step in col 4 
    diag_new( 1:ITER , 5*timept ) = momentum(1:ITER); % search momentum in col 5

    % Save files
    writematrix(diag_new, master_file,'Delimiter','tab');

return