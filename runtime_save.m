function runtime_save(runtime, testcase, E0, timept, lambda, time_start)
    
    runtime_file = [pwd '/data/runtime/runtime' testcase '_E0(' num2str(E0) ')_lambda(' num2str(lambda) ').dat'];
    if timept ~= time_start
        runtime_data = readmatrix(runtime_file);
        runtime_update = [ runtime_data ; runtime ];
    else
        runtime_update = runtime;
    end
    writematrix(runtime_update, runtime_file,'Delimiter','tab');

return