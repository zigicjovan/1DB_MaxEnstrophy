function branch_save(timept,f_ens,E0,lambda,T,testcase)
%{
Add branch evolution ( For each E0 at T = 1:end) of Enstrophy and Energy:
(1) [ Time, Max Enstrophy ]
<=>
[ T , max(E_T(phi^(n))) ]
(2) [ Time, Max Energy ]
<=>
[ T , max(K_T(phi^(n))) ]
%}

    branch_maxE0_file = [pwd '/data/enstrophy_solution/maxenstrophy'...
        testcase '_E0(' num2str(E0) ')_lambda(' num2str(lambda) ')_case3.dat'];
    branch_finalE0_file = [pwd '/data/enstrophy_solution/finalenstrophy'...
        testcase '_E0(' num2str(E0) ')_lambda(' num2str(lambda) ')_case3.dat'];

    % % Make final energy evolution file
    % % [ T , max(K_T(phi^(n))) ]
    % energy_evol_final = zeros(1, 2); % make 2 column matrix
    % energy_evol_final(1,1) = T; % time point number in col 1
    % energy_evol_final(1,2) = K(ITER); % final energy in col 2

    % Make final enstrophy evolution file
    % [ T , max(E_T(phi^(n))) , branch point , time point ]
    enstrophy_branch_final = zeros(1, 4); % make 2 column matrix
    enstrophy_branch_final(1,1) = T; % time point number in col 1
    enstrophy_branch_final(1,2) = f_ens( end ); % final enstrophy in col 2
    enstrophy_branch_final(1,3) = timept; % branch point in col 3
    enstrophy_branch_final(1,4) = length( f_ens ); % final enstrophy time point in col 4

    % Make max enstrophy evolution file
    % [ T , max(E(T)) , branch point , time point ]
    enstrophy_branch_max = zeros(1, 4); % make 2 column matrix
    enstrophy_branch_max(1,1) = T; % time point number in col 1
    [ maxens , indmax ] = max(f_ens);
    enstrophy_branch_max(1,2) = maxens; % max enstrophy in col 2
    enstrophy_branch_max(1,3) = timept; % branch point in col 3
    enstrophy_branch_max(1,4) = indmax; % max enstrophy time point in col 4
    disp(['Max Enstrophy ' num2str(enstrophy_branch_max(1,2)) ' at timepoint ' num2str(enstrophy_branch_max(1,4)) ]);

    % Save files (append for existing branch)
    writematrix(enstrophy_branch_max, branch_maxE0_file,'WriteMode','append','Delimiter','tab');
    writematrix(enstrophy_branch_final, branch_finalE0_file,'WriteMode','append','Delimiter','tab');

return