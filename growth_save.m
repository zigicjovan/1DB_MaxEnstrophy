function growth_save(lambda,E0,testcase)
    
    E0_vsEmax_file = [pwd '/data/enstrophy_solution/E0_vsEmax' testcase '_lambda(' num2str(lambda) ')_case3.dat'];

    % Make optimal IC for [E0 , T] file
    % [ E0, E_max, index of max T, index of max t , max T ]
    branch_maxE0_file = [pwd '/data/enstrophy_solution/maxenstrophy' testcase '_E0(' num2str(E0) ')_lambda(' num2str(lambda) ')_case3.dat'];
    enstrophy_branch_max = readmatrix(branch_maxE0_file);
    E0_vsEmax = zeros(1, 5); % make 4 column matrix
    E0_vsEmax(1,1) = E0; % initial enstrophy in col 1
    [ maxens , indmax ] = max(enstrophy_branch_max(:,2));
    E0_vsEmax(1,2) = maxens + E0; 
    E0_vsEmax(1,3) = enstrophy_branch_max(indmax,3); 
    E0_vsEmax(1,4) = enstrophy_branch_max(indmax,4); 
    E0_vsEmax(1,5) = enstrophy_branch_max(indmax,1); 

    % Save files
    writematrix(E0_vsEmax, E0_vsEmax_file,'WriteMode','append','Delimiter','tab'); 

return