    
Enstrophy = logspace(3,6,31);
enspts = 15;
lambda = 1.02;
casename = 'RHB_L280_K0.94b'; 
E0_vsEmax_file = [pwd '/data/enstrophy_solution/E0_vsEmax_2048_LuMax_t1_' casename '_lambda(' num2str(lambda) ').dat'];

Emax = zeros(enspts,5);

for enspt = 1:enspts

    timepts = 16;
    E0 = Enstrophy(enspt);
    Emax(enspt,1) = E0;
    branch_maxE0_file_0 = [pwd '/data/enstrophy_solution/maxenstrophy_2048_LuCont0_t1_' casename '_E0(' num2str(E0) ')_lambda(' num2str(lambda) ').dat'];
    branch_maxE0_file_1 = [pwd '/data/enstrophy_solution/maxenstrophy_2048_LuCont1_t1_' casename '_E0(' num2str(E0) ')_lambda(' num2str(lambda) ').dat'];
    branch_maxE0_file_max = [pwd '/data/enstrophy_solution/maxenstrophy_2048_LuMax_t1_' casename '_E0(' num2str(E0) ')_lambda(' num2str(lambda) ').dat'];

    enstrophy_branch_0 = readmatrix(branch_maxE0_file_0);
    enstrophy_branch_1 = readmatrix(branch_maxE0_file_1);

    enstrophy_branch_max = zeros(timepts,size(enstrophy_branch_0,2));

    for j = 1:timepts
        if enstrophy_branch_0(j,2) > enstrophy_branch_1(j,2)
            enstrophy_branch_max(j,:) = enstrophy_branch_0(j,:);
        else
            enstrophy_branch_max(j,:) = enstrophy_branch_1(j,:);
        end
    end

    maxrow = find( enstrophy_branch_max(:,2) == max(enstrophy_branch_max(:,2)) , 1 ,  'first' ) ;
    Emax(enspt,2:4) = enstrophy_branch_max(maxrow,2:4);
    Emax(enspt,5) = enstrophy_branch_max(maxrow,1);
    writematrix(enstrophy_branch_max, branch_maxE0_file_max,'Delimiter','tab');

end

writematrix(Emax, E0_vsEmax_file,'Delimiter','tab');

