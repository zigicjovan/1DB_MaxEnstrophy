function kappa_save(kappa,epsilon,timept,E0,lambda,ITER,testcase)
%{
Test Riesz representation theorem of Gateaux derivative to L_2 inner
product
%}

    kappa_E0_file = [pwd '/data/kappa/kappa_E0_' num2str(E0) '/kappa' testcase '_E0(' num2str(E0) ')' ...
        '_Timept_' num2str(timept) '_lambda(' num2str(lambda) ').dat'];
    abskappa_E0_file = [pwd '/data/kappa/kappa_E0_' num2str(E0) '/abskappa' testcase '_E0(' num2str(E0) ')' ...
        '_Timept_' num2str(timept) '_lambda(' num2str(lambda) ').dat'];

    % Make kappa file
    % [ epsilon, kappa(iter = 1), kappa(iter = 2), ... ]
    kappa_vs_eps = zeros(length(epsilon), ITER);
    kappa_vs_eps(:,ITER) = kappa';

    % Make abskappa file
    % [ epsilon, | 1 - kappa(iter = 1) |, | 1 - kappa(iter = 2) |, ... ]
    abskappa_vs_eps = zeros(length(epsilon), ITER);
    abskappa_vs_eps(:,ITER) = abs(ones(length(kappa),1) - kappa);

    if ITER > 2
        current_kappa = readmatrix(kappa_E0_file);
        current_abskappa = readmatrix(abskappa_E0_file);
        kappa_vs_eps = [ current_kappa , kappa_vs_eps(:,ITER) ];
        abskappa_vs_eps = [ current_abskappa , abskappa_vs_eps(:,ITER) ];
    else
        kappa_vs_eps(:,1) = epsilon'; 
        abskappa_vs_eps(:,1) = epsilon'; 
    end

    % Save updated files
    writematrix(kappa_vs_eps, kappa_E0_file,'Delimiter','tab');
    writematrix(abskappa_vs_eps, abskappa_E0_file,'Delimiter','tab');

return