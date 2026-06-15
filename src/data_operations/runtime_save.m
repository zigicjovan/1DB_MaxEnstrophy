function runtime_save(runtime, testcase, E0, timept, lambda)                                                        %

    runtime_file = [pwd '/data/runtime/runtime' testcase '_E0(' num2str(E0) ')_lambda(' num2str(lambda) ').dat'];   %

    % append new rows   
    try                                         %
        current = readmatrix(runtime_file);     %
        currentsize = size(current,1);          %
        
        if currentsize >= timept            %
            runtime_update = current;       %
        else
            runtime_update = NaN( timept , 2 );                 %
            runtime_update( 1:currentsize , : ) = current;      %
        end
    catch
        runtime_update = NaN(timept, 2);        % make 3 column matrix
    end
    
    % Make master runtime file
    %[ timept, runtime, cumulative ]
    runtime_update(timept,1) = timept;          % timept number in col 1
    runtime_update(timept,2) = runtime;         % timept runtime in col 2

    % Save files
    writematrix(runtime_update, runtime_file,'Delimiter','tab');    %

return