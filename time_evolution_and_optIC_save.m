function f_ens = time_evolution_and_optIC_save(Ntime,phi,timept,E0,lambda,T,N,uField,testcase,K0,stepscale)
%{
Part 1: Add time evolution ( For each E0 at T = 1, 0 mod 10, end) of Enstrophy:
(1) [ Time, Enstrophy(time) ]
<=>
[ T , E(T) ]
Part 2: Add optimal initial conditions in Fourier and Physical Space:
(1) [ wave number, abs val of fourier coeff ]
<=>
[ k , fft(phi) ]
(2) [ Grid points, Optimal IC in Physical Space ]
<=>
[ N, interpft(real(ifft(phi)) ]
%}

    terminal_E0_file = [pwd '/data/time_evolution/terminal' testcase...
        '_E0(' num2str(E0) ')_lambda(' num2str(lambda) ').dat'];

    optIC_phys_file = [pwd '/data/spectrum/spectrum_E0_' num2str(E0)...
        '/optICphys' testcase '_E0(' num2str(E0) ')_Timept_'...
        num2str(timept) '_lambda(' num2str(lambda) ').dat'];
    optIC_four_file = [pwd '/data/spectrum/spectrum_E0_' num2str(E0)...
        '/optICfour' testcase '_E0(' num2str(E0) ')_Timept_'...
        num2str(timept) '_lambda(' num2str(lambda) ').dat'];
    optIC_fourfull_file = [pwd '/data/spectrum/spectrum_E0_' num2str(E0)...
        '/optICfourfull' testcase '_E0(' num2str(E0) ')_Timept_'...
        num2str(timept) '_lambda(' num2str(lambda) ').dat'];
    deriv_file = [pwd '/data/spectrum/spectrum_E0_' num2str(E0) '/deriv_optICphys'...
        testcase '_E0(' num2str(E0) ')_Timept_'...
        num2str(timept) '_lambda(' num2str(lambda) ').dat'];
    derivfdm_file = [pwd '/data/spectrum/spectrum_E0_' num2str(E0) '/derivfdm_optICphys'...
        testcase '_E0(' num2str(E0) ')_Timept_'...
        num2str(timept) '_lambda(' num2str(lambda) ').dat'];

    % Make enstrophy time growth file
    % track spectrum for enstrophy time evolution
    eval_pts = ceil(Ntime/stepscale);
    t_evolution = linspace(0,T,Ntime);
    dt = t_evolution(2) - t_evolution(1);
    f_ens = zeros(size(t_evolution));
    K_w = N/2;
    optPhi_phys = zeros(N, eval_pts);
    du_phys = zeros(N, eval_pts); % absolute derivative of physical solution
    du_phys_fdm = zeros(N, eval_pts); 
    optPhi_four = zeros(K_w, eval_pts);
    optPhi_fourfull = zeros(N, eval_pts);
    optPhi_phys(:,1) = 1:N; % physical space
    du_phys(:,1) = optPhi_phys(:,1); 
    du_phys_fdm(:,1) = optPhi_phys(:,1); 
    optPhi_four(:,1) = 1:K_w; % wavenumbers
    optPhi_fourfull(:,1) = 1:N; % wavenumbers

    for i = 2:eval_pts 
        i_current = ceil(length(f_ens)*(i/eval_pts));
        ut = uField(i_current,:);
        f_ens(i_current) = eval_J(ut,phi,K0,N); % enstrophy value at t_i
        f_phi = adjust_optIC( phi , E0 + f_ens(i_current) , K0 , N ); % phi value at t_i  
        f_phi_x = (2*pi*1i*K0).*f_phi; % spectral derivative
        optPhi_phys(:,i) = interpft(real(ifft(f_phi)),N);
        du_phys(:,i) = abs(real(ifft(f_phi_x))); % absolute value of derivative

        du_phys_fdm(1,i) = abs( ( optPhi_phys(2,i) - 2*optPhi_phys(1,i) + optPhi_phys(N,i) ) / (dt^2) ) ;
        for j = 2 : N-1
            du_phys_fdm(j,i) = abs( ( optPhi_phys(j+1,i) - 2*optPhi_phys(j,i) + optPhi_phys(j-1,i) ) / (dt^2) ) ;
        end
        du_phys_fdm(N,i) = abs( ( optPhi_phys(1,i) - 2*optPhi_phys(N,i) + optPhi_phys(N-1,i) ) / (dt^2) ) ;

        % shift periodic solution
        if i == 2 
            shift = find( optPhi_phys(:,i) < 1e-1 & optPhi_phys(:,i) > -1e-1 , 1, 'first');
        end
        if shift > 1
            optPhi_phys(:,i) = [ optPhi_phys(shift:end,i) ; optPhi_phys(1:shift-1,i) ];
            du_phys(:,i) = [ du_phys(shift:end,i) ; du_phys(1:shift-1,i) ];
            du_phys_fdm(:,i) = [ du_phys_fdm(shift:end,i) ; du_phys_fdm(1:shift-1,i) ];
        end 
        optPhi_four(:,i) = abs(f_phi(2 : K_w + 1));
        optPhi_fourfull(:,i) = f_phi;
    end

    % Save spectrum files
    writematrix(optPhi_phys, optIC_phys_file,'Delimiter','tab');
    writematrix(optPhi_four, optIC_four_file,'Delimiter','tab');
    writematrix(optPhi_fourfull, optIC_fourfull_file,'Delimiter','tab');
    writematrix(du_phys, deriv_file,'Delimiter','tab');
    writematrix(du_phys_fdm, derivfdm_file,'Delimiter','tab');

    t_evolution( : , all(~f_ens,1) ) = []; % remove zero cols corresponding to f_ens
    t_evolution = [ 0, t_evolution ]; % start at t = 0
    f_ens( : , all(~f_ens,1) ) = []; % remove zero cols corresponding to f_ens
    f_ens = [ 0, f_ens ]; % start at ET - E0 = 0
    current = zeros(1);
    try
        current = readmatrix(terminal_E0_file);
    catch
        %
    end
    current_size = [ size(current,1) , size(f_ens,2) ];
    size_choice = max(current_size);
    enstrophy_time_final = zeros(size_choice, 2*timept); % make 2 column matrix
    try
        current = readmatrix(terminal_E0_file);
        if size(current,2) >= size(enstrophy_time_final,2)
            enstrophy_time_final = readmatrix(terminal_E0_file);
        elseif size(current,2) < size(enstrophy_time_final,2)
            enstrophy_time_final( 1:size(current,1) , 1:((2*timept)-2) ) = readmatrix(terminal_E0_file); 
        end
    catch
        %
    end
    enstrophy_time_final( 1:size(t_evolution,2) , (2*timept)-1 ) = t_evolution; % append new time window
    enstrophy_time_final( 1:size(f_ens,2) , 2*timept ) = f_ens; % enstrophy in appended column
    enstrophy_time_final( enstrophy_time_final == 0 ) = NaN; % replace zeros with NaN
    enstrophy_time_final( 1 , : ) = 0; % replace first row NaN with zeros

    writematrix(enstrophy_time_final, terminal_E0_file,'Delimiter','tab');

return