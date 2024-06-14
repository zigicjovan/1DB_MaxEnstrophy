function evaluation_save(J,timept,lambda,E0,ITER,testcase)
%{
(1) [ objective functional evaluation iteration, objective functional value ]
<=>
[ n(iter), J(iter) ]
%}
    
    % Make file
    %[ n(iter), J(iter) ]
    diagnostics = zeros(ITER, 2); % make 2 column matrix
    diagnostics(:,1) = 1:ITER; % iteration number in col 1
    diagnostics(:,2) = J(1:ITER); % objective functional value in col 2

    master_file = [pwd '/data/diagnostics/objective' testcase '_E0(' num2str(E0) ')_lambda(' num2str(lambda) ').dat'];

    % append new columns
    try
        current = readmatrix(master_file);

        current_size = [ size(current,1) , size(diagnostics,1) ];
        size_choice = max(current_size);
        diag_new = NaN(size_choice, 2*timept); % make 2 column matrix
        
        if size(current,2) >= size(diag_new,2)
            diag_new = readmatrix(master_file);
        elseif size(current,2) < size(diag_new,2)
            diag_new( 1:size(current,1) , 1:((2*timept)-2) ) = readmatrix(master_file); 
        end

        % append new values
        diag_new( 1:size(diagnostics,1) , (2*timept)-1 ) = diagnostics(:,1); 
        diag_new( 1:size(diagnostics,1) , 2*timept ) = diagnostics(:,2); 
        diagnostics = diag_new;
    catch
        %
    end

    % Save files
    writematrix(diagnostics, master_file,'Delimiter','tab');

return