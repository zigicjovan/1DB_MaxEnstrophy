function evaluation_save(J,timept,lambda,E0,ITER,testcase)
%{
(1) [ objective functional evaluation iteration, objective functional value ]
<=>
[ n(iter), J(iter) ]
%}

    master_file = [pwd '/data/diagnostics/objective' testcase '_E0(' num2str(E0) ')_lambda(' num2str(lambda) ').dat'];

    % append new columns
    try
        current = readmatrix(master_file);
        current_size = [ size(current,1) , ITER ];
        size_choice = max(current_size);
        diag_new = NaN(size_choice, 2*timept); % make 2 column matrix
        diag_new( 1:size(current,1) , 1:size(current,2) ) = current; 
    catch
        diag_new = NaN(ITER, 2*timept); % make 2 column matrix
    end

    % append new values
    diag_new( 1:ITER , (2*timept)-1 ) = 1:ITER; % iteration number in col 1
    diag_new( 1:ITER , 2*timept ) = J(1:ITER); % objective functional value in col 2

    % Save files
    writematrix(diag_new, master_file,'Delimiter','tab');

return