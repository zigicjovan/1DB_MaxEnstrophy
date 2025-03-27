function max_save(timept,lambda,E0,testcase_max,testcase_old,testcase_other,timeend)

    %%% diagnostics 
    master_file_old = [pwd '/data/diagnostics/diag' testcase_old '_E0(' num2str(E0) ')_lambda(' num2str(lambda) ').dat'];

    master_file_other = [pwd '/data/diagnostics/diag' testcase_other '_E0(' num2str(E0) ')_lambda(' num2str(lambda) ').dat'];

    master_file_new = [pwd '/data/diagnostics/diag' testcase_max '_E0(' num2str(E0) ')_lambda(' num2str(lambda) ').dat'];

    objective_file_old = [pwd '/data/diagnostics/objective' testcase_old '_E0(' num2str(E0) ')_lambda(' num2str(lambda) ').dat'];

    objective_file_other = [pwd '/data/diagnostics/objective' testcase_other '_E0(' num2str(E0) ')_lambda(' num2str(lambda) ').dat'];

    objective_file_new = [pwd '/data/diagnostics/objective' testcase_max '_E0(' num2str(E0) ')_lambda(' num2str(lambda) ').dat'];

    %%% runtime 
    runtime_file_old = [pwd '/data/runtime/runtime' testcase_old ...
        '_E0(' num2str(E0) ')_lambda(' num2str(lambda) ').dat'];

    runtime_file_other = [pwd '/data/runtime/runtime' testcase_other ...
        '_E0(' num2str(E0) ')_lambda(' num2str(lambda) ').dat'];

    runtime_file_new = [pwd '/data/runtime/runtime' testcase_max ...
        '_E0(' num2str(E0) ')_lambda(' num2str(lambda) ').dat'];

    %%% time evol and opt IC
    terminal_E0_file_old = [pwd '/data/time_evolution/terminal' testcase_old...
        '_E0(' num2str(E0) ')_lambda(' num2str(lambda) ').dat'];
    optIC_phys_file_old = [pwd '/data/spectrum/temp/optICphys'...
        testcase_old '_E0(' num2str(E0) ')_Timept_'...
        num2str(timept) '_lambda(' num2str(lambda) ').dat'];
    optIC_four_file_old = [pwd '/data/spectrum/temp/optICfour'...
        testcase_old '_E0(' num2str(E0) ')_Timept_'...
        num2str(timept) '_lambda(' num2str(lambda) ').dat'];
    %{
    optIC_fourfull_file_old = [pwd '/data/spectrum/temp/optICfourfull'...
        testcase_old '_E0(' num2str(E0) ')_Timept_'...
        num2str(timept) '_lambda(' num2str(lambda) ').dat'];
    deriv_file_old = [pwd '/data/spectrum/temp/deriv_optICphys'...
        testcase_old '_E0(' num2str(E0) ')_Timept_'...
        num2str(timept) '_lambda(' num2str(lambda) ').dat'];
    derivfdm_file_old = [pwd '/data/spectrum/temp/derivfdm_optICphys'...
        testcase_old '_E0(' num2str(E0) ')_Timept_'...
        num2str(timept) '_lambda(' num2str(lambda) ').dat'];
    %}

    terminal_E0_file_other = [pwd '/data/time_evolution/terminal' testcase_other...
        '_E0(' num2str(E0) ')_lambda(' num2str(lambda) ').dat'];
    optIC_phys_file_other = [pwd '/data/spectrum/temp/optICphys'...
        testcase_other '_E0(' num2str(E0) ')_Timept_'...
        num2str(timept) '_lambda(' num2str(lambda) ').dat'];
    optIC_four_file_other = [pwd '/data/spectrum/temp/optICfour'...
        testcase_other '_E0(' num2str(E0) ')_Timept_'...
        num2str(timept) '_lambda(' num2str(lambda) ').dat'];
    %{
    optIC_fourfull_file_other = [pwd '/data/spectrum/temp/optICfourfull'...
        testcase_other '_E0(' num2str(E0) ')_Timept_'...
        num2str(timept) '_lambda(' num2str(lambda) ').dat'];
    deriv_file_other = [pwd '/data/spectrum/temp/deriv_optICphys'...
        testcase_other '_E0(' num2str(E0) ')_Timept_'...
        num2str(timept) '_lambda(' num2str(lambda) ').dat'];
    derivfdm_file_other = [pwd '/data/spectrum/temp/derivfdm_optICphys'...
        testcase_other '_E0(' num2str(E0) ')_Timept_'...
        num2str(timept) '_lambda(' num2str(lambda) ').dat'];
    %}

    terminal_E0_file_new = [pwd '/data/time_evolution/terminal' testcase_max...
        '_E0(' num2str(E0) ')_lambda(' num2str(lambda) ').dat'];
    optIC_phys_file_new = [pwd '/data/spectrum/spectrum_E0_' num2str(E0)...
        '/optICphys' testcase_max '_E0(' num2str(E0) ')_Timept_'...
        num2str(timept) '_lambda(' num2str(lambda) ').dat'];
    optIC_four_file_new = [pwd '/data/spectrum/spectrum_E0_' num2str(E0)...
        '/optICfour' testcase_max '_E0(' num2str(E0) ')_Timept_'...
        num2str(timept) '_lambda(' num2str(lambda) ').dat'];
    %{
    optIC_fourfull_file_new = [pwd '/data/spectrum/spectrum_E0_' num2str(E0)...
        '/optICfourfull' testcase_max '_E0(' num2str(E0) ')_Timept_'...
        num2str(timept) '_lambda(' num2str(lambda) ').dat'];
    deriv_file_new = [pwd '/data/spectrum/spectrum_E0_' num2str(E0) '/deriv_optICphys'...
        testcase_max '_E0(' num2str(E0) ')_Timept_'...
        num2str(timept) '_lambda(' num2str(lambda) ').dat'];
    derivfdm_file_new = [pwd '/data/spectrum/spectrum_E0_' num2str(E0) '/derivfdm_optICphys'...
        testcase_max '_E0(' num2str(E0) ')_Timept_'...
        num2str(timept) '_lambda(' num2str(lambda) ').dat'];
    %}

    %%% branch
    branch_maxE0_file_old = [pwd '/data/enstrophy_solution/maxenstrophy'...
        testcase_old '_E0(' num2str(E0) ')_lambda(' num2str(lambda) ').dat'];
    branch_maxE0_file_other = [pwd '/data/enstrophy_solution/maxenstrophy'...
        testcase_other '_E0(' num2str(E0) ')_lambda(' num2str(lambda) ').dat'];
    branch_maxE0_file_new = [pwd '/data/enstrophy_solution/maxenstrophy'...
        testcase_max '_E0(' num2str(E0) ')_lambda(' num2str(lambda) ').dat'];

    %{
    branch_finalE0_file_old = [pwd '/data/enstrophy_solution/finalenstrophy'...
        testcase_old '_E0(' num2str(E0) ')_lambda(' num2str(lambda) ').dat'];
    branch_finalE0_file_other = [pwd '/data/enstrophy_solution/finalenstrophy'...
        testcase_other '_E0(' num2str(E0) ')_lambda(' num2str(lambda) ').dat'];
    branch_finalE0_file_new = [pwd '/data/enstrophy_solution/finalenstrophy'...
        testcase_max '_E0(' num2str(E0) ')_lambda(' num2str(lambda) ').dat'];
    %}

    % Load files
    diagnostics = readmatrix(master_file_old);
    objective = readmatrix(objective_file_old);

    runtime = readmatrix(runtime_file_old);
    runtime2 = readmatrix(runtime_file_other);
    currentrow = find( runtime(:,1) == timept , 1 ,  'last' ) ;
    try
        runtime_partimept = max( runtime(currentrow,2) , runtime2(currentrow,2) );
        runtime_seqtimept = runtime(currentrow,2) + runtime2(currentrow,2) ;
        runtime_new = readmatrix(runtime_file_new);
        runtime_new( currentrow , 1:3 ) = runtime( currentrow , 1:3 ); 
        runtime = runtime_new;
        runtime(currentrow,4) = runtime_partimept; % parallel time
        runtime(currentrow,6) = runtime_seqtimept; % sequential time
        if isnan(runtime(currentrow-1,3))
            runtime(currentrow,3) = runtime(currentrow,2); 
            runtime(currentrow,5) = runtime(currentrow,4);
            runtime(currentrow,7) = runtime(currentrow,6);
        else
            runtime(currentrow,3) = runtime(currentrow,2) + runtime(currentrow-1,3); % cumulative max
            runtime(currentrow,5) = runtime(currentrow,4) + runtime(currentrow-1,5); % cumulative parallel time
            runtime(currentrow,7) = runtime(currentrow,6) + runtime(currentrow-1,7); % cumulative sequential time
        end
    catch
        runtime_partimept = max( runtime(currentrow,2) , runtime2(currentrow,2) );
        runtime_seqtimept = runtime(currentrow,2) + runtime2(currentrow,2) ;
        runtime_new = NaN(currentrow,7);
        runtime_new(currentrow,1:3) = runtime( currentrow , 1:3 );
        runtime_new(currentrow,4) = runtime_partimept; % parallel time
        runtime_new(currentrow,5) = runtime_partimept; % cumulative parallel time
        runtime_new(currentrow,6) = runtime_seqtimept; % sequential time
        runtime_new(currentrow,7) = runtime_seqtimept; % cumulative sequential time
        runtime = runtime_new;
    end

    enstrophy_time_final = readmatrix(terminal_E0_file_old);
    optPhi_phys = readmatrix(optIC_phys_file_old);
    optPhi_four = readmatrix(optIC_four_file_old);
    %{
    optPhi_fourfull = readmatrix(optIC_fourfull_file_old);
    du_phys = readmatrix(deriv_file_old);
    du_phys_fdm = readmatrix(derivfdm_file_old);
    %}

    enstrophy_branch_max = readmatrix(branch_maxE0_file_old);
    % enstrophy_branch_final = readmatrix(branch_finalE0_file_old);
    currentrow = find( enstrophy_branch_max(:,3) == timept , 1 ,  'last' ) ;
    try
        enstrophy_branch_max_new = readmatrix(branch_maxE0_file_new);
        enstrophy_branch_max_new( currentrow, : ) = enstrophy_branch_max( currentrow , : );
        enstrophy_branch_max = enstrophy_branch_max_new;
        %{
        enstrophy_branch_final_new = readmatrix(branch_finalE0_file_new);
        enstrophy_branch_final_new( timept , : ) = enstrophy_branch_final( timept , : );
        enstrophy_branch_final = enstrophy_branch_final_new; 
        %}
    catch
        enstrophy_branch_max_new = NaN( currentrow , size(enstrophy_branch_max,2) );
        enstrophy_branch_max_new( currentrow , : ) = enstrophy_branch_max( currentrow , : );
        enstrophy_branch_max = enstrophy_branch_max_new;
        %{
        enstrophy_branch_final_new = zeros( timept , size(enstrophy_branch_final,2) );
        enstrophy_branch_final_new ( timept , : )= enstrophy_branch_final( 1 , : );
        enstrophy_branch_final = enstrophy_branch_final_new; 
        %}
    end

    try
        enstrophy_time_final_new = readmatrix(terminal_E0_file_new);

        r_new = size(enstrophy_time_final_new,1);
        r = size(enstrophy_time_final,1);
        c_new = size(enstrophy_time_final_new,2);
        c = size(enstrophy_time_final,2);

        % append to evolution files
        if r_new < r
            enstrophy_time_final( 1:r_new , 1:c_new ) = enstrophy_time_final_new( 1:r_new , 1:c_new );
            enstrophy_time_final( r_new+1:r , 1:c_new ) = NaN( r-r_new , c_new );
        else
            enstrophy_time_final_x = NaN(r_new,c);
            enstrophy_time_final_x( 1:r_new , 1:c_new ) = enstrophy_time_final_new( 1:r_new , 1:c_new );
            enstrophy_time_final_x( 1:r , (c_new+1):c ) = enstrophy_time_final( 1:r , (c_new+1):c );
            enstrophy_time_final = enstrophy_time_final_x;
        end

    catch
        %
    end

    try
        diagnostics_new = readmatrix(master_file_new);

        r2_new = size(diagnostics_new,1);
        r2 = size(diagnostics,1);
        c_new = size(diagnostics_new,2);
        c = size(diagnostics,2);

        % append 
        if r2_new < r2
            diagnostics( 1:r2_new , 1:c_new ) = diagnostics_new( 1:r2_new , 1:c_new );
            diagnostics( r2_new+1:r2 , 1:c_new ) = NaN( r2-r2_new , c_new );
        else
            diagnostics_x = NaN(r2_new,c);
            diagnostics_x( 1:r2_new , 1:c_new ) = diagnostics_new( 1:r2_new , 1:c_new );
            diagnostics_x( 1:r2 , (c_new+1):c ) = diagnostics( 1:r2 , (c_new+1):c );
            diagnostics = diagnostics_x;
        end

    catch
        %
    end

    try
        objective_new = readmatrix(objective_file_new);

        r2_new = size(objective_new,1);
        r2 = size(objective,1);

        % append 
        if r2_new < r2
            objective( 1:r2_new , 1:(2*timept)-2 ) = objective_new( 1:r2_new , 1:(2*timept)-2 );
            objective( r2_new+1:r2 , 1:(2*timept)-2 ) = NaN( r2-r2_new , (2*timept)-2 );
        else
            objective_x = NaN(r2_new,2*timept);
            objective_x( 1:r2_new , 1:(2*timept)-2 ) = objective_new( 1:r2_new , 1:(2*timept)-2 );
            objective_x( 1:r2 , (2*timept)-1:2*timept ) = objective( 1:r2 , (2*timept)-1:2*timept );
            objective = objective_x;
        end

    catch
        %
    end

    % Save files
    writematrix(diagnostics, master_file_new,'Delimiter','tab');
    writematrix(objective, objective_file_new,'Delimiter','tab');
    writematrix(runtime, runtime_file_new,'Delimiter','tab');

    writematrix(enstrophy_time_final, terminal_E0_file_new,'Delimiter','tab');
    writematrix(optPhi_phys, optIC_phys_file_new,'Delimiter','tab');
    writematrix(optPhi_four, optIC_four_file_new,'Delimiter','tab');
    %{
    writematrix(optPhi_fourfull, optIC_fourfull_file_new,'Delimiter','tab');
    writematrix(du_phys, deriv_file_new,'Delimiter','tab');
    writematrix(du_phys_fdm, derivfdm_file_new,'Delimiter','tab');
    %}

    writematrix(enstrophy_branch_max, branch_maxE0_file_new,'Delimiter','tab');
    % writematrix(enstrophy_branch_final, branch_finalE0_file_new,'Delimiter','tab');

    % Delete timept files (not terminal or branch) and delete terminal or branch at last timept
    delete(optIC_phys_file_old);
    delete(optIC_four_file_old);
    %{
    delete(optIC_fourfull_file_old);
    delete(deriv_file_old);
    delete(derivfdm_file_old);
    %}
    
    delete(optIC_phys_file_other);
    delete(optIC_four_file_other);
    %{
    delete(optIC_fourfull_file_other);
    delete(deriv_file_other);
    delete(derivfdm_file_other);
    %}

    if timept == timeend % && timept == 999
        delete(master_file_old);
        delete(objective_file_old);
        delete(runtime_file_old);  
        delete(terminal_E0_file_old);
        delete(branch_maxE0_file_old);
        % delete(branch_finalE0_file_old);
        
        delete(master_file_other);
        delete(objective_file_other);
        delete(runtime_file_other);
        delete(terminal_E0_file_other);
        delete(branch_maxE0_file_other);
        % delete(branch_finalE0_file_other);
    end

return