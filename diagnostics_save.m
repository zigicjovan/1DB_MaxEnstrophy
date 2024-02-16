function diagnostics_save(J,K,step,momentum,timept,lambda,E0,ITER,testcase)
%{
VERBOSE 1: Add master diagnostics file and optimal initial conditions in physical space:
(1) [ Iterations, Enstrophy, Energy, Search Step Size, Search Momentum ]
<=>
[ n(iter), E_T(phi^(n)), K_T(phi^(n)), tau_n, beta_n ]
%}

    master_file = [pwd '/data/diagnostics/diagnostics_E0_' num2str(E0)...
        '/diag' testcase '_E0(' num2str(E0) ')_Timept_' num2str(timept)...
        '_lambda(' num2str(lambda) ').dat'];

    % Make master diagnostics file
    %[ n(iter), E_T(phi^(n)), K_T(phi^(n)), tau_n, beta_n ]
    diagnostics = zeros(ITER, 5); % make 5 column matrix
    diagnostics(:,1) = 1:ITER; % iteration number in col 1
    diagnostics(:,2) = J(1:ITER); % enstrophy in col 2
    diagnostics(:,3) = K(1:ITER); % energy in col 3
    diagnostics(:,4) = step(1:ITER); % search step in col 4
    diagnostics(:,5) = momentum(1:ITER); % search momentum in col 5

    % Save files
    writematrix(diagnostics, master_file,'Delimiter','tab');

return