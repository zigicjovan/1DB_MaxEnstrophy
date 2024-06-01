function runtime_save(runtime, testcase, E0, timept, lambda)

    % Make master runtime file
    %[ timept, runtime ]
    runtime_update = zeros(1, 2); % make 5 column matrix
    runtime_update(:,1) = timept; % iteration number in col 1
    runtime_update(:,2) = runtime; % objective functional value in col 2

    runtime_file = [pwd '/data/runtime/runtime' testcase '_E0(' num2str(E0) ')_lambda(' num2str(lambda) ').dat'];

    % append new rows   
    try
        current = readmatrix(runtime_file);
        currentsize = size(current,1);
        
        if currentsize >= timept
            runtime_new = current;
        else
            runtime_new = NaN( currentsize + 1 , 2 );
            runtime_new( 1:currentsize , : ) = current; 
        end

        % append new values
        runtime_new( timept , : ) = runtime_update;
        runtime_update = runtime_new;
    catch
        %
    end

    % Save files
    writematrix(runtime_update, runtime_file,'Delimiter','tab');

return