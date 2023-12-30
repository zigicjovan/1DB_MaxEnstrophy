
Enstrophy = logspace(3,6,31);
EnstPoints = 21;
ens_start = 11;

TimePoints = 31;

testcase = '_N4096';
lambda = 0.5;

for enspt = ens_start:EnstPoints

    E0 = Enstrophy(enspt);
    prefactor = 2;
    T_ens_UB = prefactor*(1/sqrt(E0)); % Adjust for decreasing T_max    
    T_ens_LB = T_ens_UB/TimePoints; % Adjust for decreasing T_max  
    TimeWindow = linspace(T_ens_LB,T_ens_UB,TimePoints); 

    terminal_E0_file = [pwd '/data/time_evolution/terminal' testcase...
        '_E0(' num2str(E0) ')_lambda(' num2str(lambda) ')_case3.dat'];
    branch_maxE0_file = [pwd '/data/enstrophy_solution/maxenstrophy'...
        testcase '_E0(' num2str(E0) ')_lambda(' num2str(lambda) ')_case3_2.dat'];
    E0_vsEmax_file = [pwd '/data/enstrophy_solution/E0_vsEmax' testcase '_lambda(' num2str(lambda) ')_case3_2.dat'];

    for timept = 1:TimePoints
                
        T = TimeWindow(timept);
        current = readmatrix(terminal_E0_file);
        f_ens = current( : , 2*timept );

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

    end

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
end