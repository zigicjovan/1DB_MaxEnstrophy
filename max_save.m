function max_save(timept,lambda,E0,testcase_max,testcase_old)

    % diagnostics 
    master_file_old = [pwd '/data/diagnostics/diagnostics_E0_' num2str(E0)...
        '/diag' testcase_old '_E0(' num2str(E0) ')_Timept_' num2str(timept)...
        '_lambda(' num2str(lambda) ').dat'];

    master_file_new = [pwd '/data/diagnostics/diagnostics_E0_' num2str(E0)...
        '/diag' testcase_max '_E0(' num2str(E0) ')_Timept_' num2str(timept)...
        '_lambda(' num2str(lambda) ').dat'];

    % time evol and opt IC
    terminal_E0_file_old = [pwd '/data/time_evolution/terminal' testcase_old...
        '_E0(' num2str(E0) ')_lambda(' num2str(lambda) ').dat'];
    optIC_phys_file_old = [pwd '/data/spectrum/temp/optICphys'...
        testcase_old '_E0(' num2str(E0) ')_Timept_'...
        num2str(timept) '_lambda(' num2str(lambda) ').dat'];
    optIC_four_file_old = [pwd '/data/spectrum/temp/optICfour'...
        testcase_old '_E0(' num2str(E0) ')_Timept_'...
        num2str(timept) '_lambda(' num2str(lambda) ').dat'];
    optIC_fourfull_file_old = [pwd '/data/spectrum/temp/optICfourfull'...
        testcase_old '_E0(' num2str(E0) ')_Timept_'...
        num2str(timept) '_lambda(' num2str(lambda) ').dat'];
    deriv_file_old = [pwd '/data/spectrum/temp/deriv_optICphys'...
        testcase_old '_E0(' num2str(E0) ')_Timept_'...
        num2str(timept) '_lambda(' num2str(lambda) ').dat'];
    derivfdm_file_old = [pwd '/data/spectrum/temp/derivfdm_optICphys'...
        testcase_old '_E0(' num2str(E0) ')_Timept_'...
        num2str(timept) '_lambda(' num2str(lambda) ').dat'];

    terminal_E0_file_new = [pwd '/data/time_evolution/terminal' testcase_max...
        '_E0(' num2str(E0) ')_lambda(' num2str(lambda) ').dat'];
    optIC_phys_file_new = [pwd '/data/spectrum/spectrum_E0_' num2str(E0)...
        '/optICphys' testcase_max '_E0(' num2str(E0) ')_Timept_'...
        num2str(timept) '_lambda(' num2str(lambda) ').dat'];
    optIC_four_file_new = [pwd '/data/spectrum/spectrum_E0_' num2str(E0)...
        '/optICfour' testcase_max '_E0(' num2str(E0) ')_Timept_'...
        num2str(timept) '_lambda(' num2str(lambda) ').dat'];
    optIC_fourfull_file_new = [pwd '/data/spectrum/spectrum_E0_' num2str(E0)...
        '/optICfourfull' testcase_max '_E0(' num2str(E0) ')_Timept_'...
        num2str(timept) '_lambda(' num2str(lambda) ').dat'];
    deriv_file_new = [pwd '/data/spectrum/spectrum_E0_' num2str(E0) '/deriv_optICphys'...
        testcase_max '_E0(' num2str(E0) ')_Timept_'...
        num2str(timept) '_lambda(' num2str(lambda) ').dat'];
    derivfdm_file_new = [pwd '/data/spectrum/spectrum_E0_' num2str(E0) '/derivfdm_optICphys'...
        testcase_max '_E0(' num2str(E0) ')_Timept_'...
        num2str(timept) '_lambda(' num2str(lambda) ').dat'];

    % branch
    branch_maxE0_file_old = [pwd '/data/enstrophy_solution/maxenstrophy'...
        testcase_old '_E0(' num2str(E0) ')_lambda(' num2str(lambda) ').dat'];
    branch_finalE0_file_old = [pwd '/data/enstrophy_solution/finalenstrophy'...
        testcase_old '_E0(' num2str(E0) ')_lambda(' num2str(lambda) ').dat'];

    branch_maxE0_file_new = [pwd '/data/enstrophy_solution/maxenstrophy'...
        testcase_max '_E0(' num2str(E0) ')_lambda(' num2str(lambda) ').dat'];
    branch_finalE0_file_new = [pwd '/data/enstrophy_solution/finalenstrophy'...
        testcase_max '_E0(' num2str(E0) ')_lambda(' num2str(lambda) ').dat'];

    % Load files
    diagnostics = readmatrix(master_file_old);

    enstrophy_time_final = readmatrix(terminal_E0_file_old);
    optPhi_phys = readmatrix(optIC_phys_file_old);
    optPhi_four = readmatrix(optIC_four_file_old);
    optPhi_fourfull = readmatrix(optIC_fourfull_file_old);
    du_phys = readmatrix(deriv_file_old);
    du_phys_fdm = readmatrix(derivfdm_file_old);

    enstrophy_branch_max = readmatrix(branch_maxE0_file_old);
    enstrophy_branch_final = readmatrix(branch_finalE0_file_old);

    if timept > 1
        enstrophy_time_final_new = readmatrix(terminal_E0_file_new);
        enstrophy_branch_max_new = readmatrix(branch_maxE0_file_new);
        enstrophy_branch_final_new = readmatrix(branch_finalE0_file_new);
        enstrophy_branch_max(1:timept-1,:) = enstrophy_branch_max_new(1:timept-1,:);
        enstrophy_branch_final(1:timept-1,:) = enstrophy_branch_final_new(1:timept-1,:);

        r_new = size(enstrophy_time_final_new,1);
        r = size(enstrophy_time_final,1);

        % append to evolution and branch files
        if r_new < r
            enstrophy_time_final( 1:r_new , 1:(2*timept)-2 ) = enstrophy_time_final_new( 1:r_new , 1:(2*timept)-2 );
            enstrophy_time_final( r_new+1:r , 1:(2*timept)-2 ) = NaN( r-r_new , (2*timept)-2 );
        else
            enstrophy_time_final_x = NaN(r_new,2*timept);
            enstrophy_time_final_x( 1:r_new , 1:(2*timept)-2 ) = enstrophy_time_final_new( 1:r_new , 1:(2*timept)-2 );
            enstrophy_time_final_x( 1:r , (2*timept)-1:2*timept ) = enstrophy_time_final( 1:r , (2*timept)-1:2*timept );
            enstrophy_time_final = enstrophy_time_final_x;
        end
        %
    end

    % Save files
    writematrix(diagnostics, master_file_new,'Delimiter','tab');

    writematrix(enstrophy_time_final, terminal_E0_file_new,'Delimiter','tab');
    writematrix(optPhi_phys, optIC_phys_file_new,'Delimiter','tab');
    writematrix(optPhi_four, optIC_four_file_new,'Delimiter','tab');
    writematrix(optPhi_fourfull, optIC_fourfull_file_new,'Delimiter','tab');
    writematrix(du_phys, deriv_file_new,'Delimiter','tab');
    writematrix(du_phys_fdm, derivfdm_file_new,'Delimiter','tab');

    writematrix(enstrophy_branch_max, branch_maxE0_file_new,'Delimiter','tab');
    writematrix(enstrophy_branch_final, branch_finalE0_file_new,'Delimiter','tab');

    % Delete other files in spectrum folder

return