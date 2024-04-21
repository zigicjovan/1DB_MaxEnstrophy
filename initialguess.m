function [ solution , T ] = initialguess(type,x,E0,a1,x_ref,timept,tptest,testcase,lambda,ss_shift)
    
    switch type
        case 'sine'
            solution = sqrt(E0)/pi*sin(2*pi*x - pi);
            T = 0;
        case 'sin2'
            solution = sqrt(E0)/(2*pi)*sin(4*pi*x - pi);
            T = 0;
        case '2 modes'
            solution = a1*sin(2*pi*x-pi) + sqrt( E0-(a1*pi)^2 )/(2*pi)*sin(4*pi*x-pi/7);
            T = 0;
        case 'exact'
            initialdata_file = [pwd '/optIC/LuLu_u_E010_nu3_N2048.dat'];
            aux = readmatrix(initialdata_file);
            if length(x) == 2048
                solution = (aux).';
            else
                solution = interp1(x_ref, aux , x, 'spline');
            end
            solution = fft(solution);
            T = 0;
        case 'optIC'
            try
                phi_guess_file = [pwd '/optIC/phi' testcase '_E0(' num2str(E0) ').dat'];
                phi_guess = readmatrix(phi_guess_file);
                solution = phi_guess(timept,:);
                solution = real(ifft(solution));
                length_ref = length(solution);
                if length(x) ~= length_ref
                    solution = interp1(x_ref, solution , x, 'spline');
                end
                solution = fft(solution);
            catch
                disp('No optIC file found, using exact solution');
                phi = initialguess('exact', x, E0, 0, x_ref,timept,tptest,testcase,lambda);
                solution = fft(phi);
            end
            T = 0;
        case 'slingshot'

            % load solution files 
            %{
            % Level 1
            optIC_phys_file = [pwd '/data/spectrum/spectrum_E0_' num2str(E0*1e0) '/optICphys_' ...
                num2str(length(x_ref)*1) '_LuFour_long_E0(' num2str(E0*1e0) ')_Timept_' num2str(tptest) '_lambda(0.5).dat'];
            % Level 2
            optIC_phys_file = [pwd '/data/spectrum/spectrum_E0_' num2str(E0*1e0) '/optICphys_' ...
                num2str(length(x_ref)*1) '_slingshotTlong_long3_E0(' num2str(E0*1e0) ')_Timept_' num2str(tptest) '_lambda(0.5).dat'];
            % Level 3
            optIC_phys_file = [pwd '/data/spectrum/spectrum_E0_' num2str(E0*1e0) '/optICphys_' ...
                num2str(length(x_ref)*1) '_slingshot2Tlong2_long10_E0(' num2str(E0*1e0) ')_Timept_' num2str(tptest) '_lambda(0.5).dat'];
            % Level 4
            optIC_phys_file = [pwd '/data/spectrum/spectrum_E0_' num2str(E0*1e0) '/optICphys_' ...
                num2str(length(x_ref)*1) '_slingshot3Tfull_long29_E0(' num2str(E0*1e0) ')_Timept_' num2str(tptest) '_lambda(0.5).dat'];
            % Level 5
            optIC_phys_file = [pwd '/data/spectrum/spectrum_E0_' num2str(E0*1e0) '/optICphys_' ...
                num2str(length(x_ref)*1) '_slingshot4Tlong_full19_E0(' num2str(E0*1e0) ')_Timept_' num2str(tptest) '_lambda(0.5).dat'];
            %}

            % find indices of abs(min) and max pts in enstrophy evolution
            term0_file = [pwd '/data/time_evolution/terminal_2048_LuMax_t1_fine50_E0(' num2str(E0) ')_lambda(' num2str(lambda) ').dat'];  
            term0 = readmatrix(term0_file);
            [ ~ , t0in ] = min( abs( term0(3:end,2*tptest) ) );
            t0in = t0in+2;
            [ ~ , t0ix ] = max(term0(:,2*tptest));
            T = term0(t0ix,(2*tptest)-1) - term0(t0in,(2*tptest)-1); % T = T(Emax) - T(root)

            % set phi IC
            optIC_four_file = [pwd '/data/spectrum/spectrum_E0_' num2str(E0) '/optICfourfull_2048_LuMax_t1_fine50_E0(' ...
                num2str(E0) ')_Timept_' num2str(tptest) '_lambda(' num2str(lambda) ').dat'];
            optIC_four = readmatrix(optIC_four_file);  
            solution = optIC_four( : , t0in );
            solution = solution';

            % Handle different resolutions
            if length(x) ~= length(x_ref)
                solution = real(ifft(solution));
                solution = interp1(x_ref, solution , x, 'spline');
                solution = fft(solution);
            end 
        case 'ss_shift'

            shift = ss_shift;

            % find indices of abs(min) and max pts in enstrophy evolution
            term0_file = [pwd '/data/time_evolution/terminal_2048_LuMax_t1_fine50_E0(' num2str(E0) ')_lambda(' num2str(lambda) ').dat'];  
            term0 = readmatrix(term0_file);
            [ ~ , t0ind ] = min( abs( term0(3:end,2*tptest) ) );
            t0in = t0ind + shift;
            t0in_root = t0ind + 2;
            [ ~ , t0ix ] = max(term0(:,2*tptest));
            T = term0(t0ix,(2*tptest)-1) - term0(t0in_root,(2*tptest)-1); % T = T(Emax) - T(root)

            % set phi IC
            optIC_four_file = [pwd '/data/spectrum/spectrum_E0_' num2str(E0) '/optICfourfull_2048_LuMax_t1_fine50_E0(' ...
                num2str(E0) ')_Timept_' num2str(tptest) '_lambda(' num2str(lambda) ').dat'];
            optIC_four = readmatrix(optIC_four_file);  
            solution = optIC_four( : , t0in );
            solution = solution';

            % Handle different resolutions
            if length(x) ~= length(x_ref)
                solution = real(ifft(solution));
                solution = interp1(x_ref, solution , x, 'spline');
                solution = fft(solution);
            end 
        case 'ss_root'

            shift = ss_shift;

            % find indices of abs(min) and max pts in enstrophy evolution
            term0_file = [pwd '/data/time_evolution/terminal_2048_LuMax_t1_fine50_E0(' num2str(E0) ')_lambda(' num2str(lambda) ').dat'];  
            term0 = readmatrix(term0_file);
            [ ~ , t0ind ] = min( term0(:,2*tptest) );
            t0in = t0ind + shift;
            t0in_root = t0in;
            [ ~ , t0ix ] = max(term0(:,2*tptest));
            T = term0(t0ix,(2*tptest)-1) - term0(t0in_root,(2*tptest)-1); % T = T(Emax) - T(root)

            % set phi IC
            optIC_four_file = [pwd '/data/spectrum/spectrum_E0_' num2str(E0) '/optICfourfull_2048_LuMax_t1_fine50_E0(' ...
                num2str(E0) ')_Timept_' num2str(tptest) '_lambda(' num2str(lambda) ').dat'];
            optIC_four = readmatrix(optIC_four_file);  
            solution = optIC_four( : , t0in );
            solution = solution';

            % Handle different resolutions
            if length(x) ~= length(x_ref)
                solution = real(ifft(solution));
                solution = interp1(x_ref, solution , x, 'spline');
                solution = fft(solution);
            end 
        case 'ss_start'

            shift = ss_shift;
            
            % find indices of abs(min) and max pts in enstrophy evolution
            term0_file = [pwd '/data/time_evolution/terminal_2048_LuMax_t1_fine50_E0(' num2str(E0) ')_lambda(' num2str(lambda) ').dat'];  
            term0 = readmatrix(term0_file);
            t0ind = 1;
            t0in = t0ind + shift;
            t0in_root = t0in;
            [ ~ , t0ix ] = max(term0(:,2*tptest));
            T = term0(t0ix,(2*tptest)-1) - term0(t0in_root,(2*tptest)-1); % T = T(Emax) - T(root)

            % set phi IC
            optIC_four_file = [pwd '/data/spectrum/spectrum_E0_' num2str(E0) '/optICfourfull_2048_LuMax_t1_fine50_E0(' ...
                num2str(E0) ')_Timept_' num2str(tptest) '_lambda(' num2str(lambda) ').dat'];
            optIC_four = readmatrix(optIC_four_file);  
            solution = optIC_four( : , t0in );
            solution = solution';

            % Handle different resolutions
            if length(x) ~= length(x_ref)
                solution = real(ifft(solution));
                solution = interp1(x_ref, solution , x, 'spline');
                solution = fft(solution);
            end 
        case 'slingshot1'

            %
            tptest = 3;
            % find indices of abs(min) and max pts in enstrophy evolution
            term0_file = [pwd '/data/time_evolution/terminal_2048_LuFour_long_E0(1000)_lambda(' num2str(lambda) ').dat'];  
            term0 = readmatrix(term0_file);
            [ ~ , t0in ] = min( abs( term0(3:end,2*tptest) ) );
            t0in = t0in+30;
            [ ~ , t0ix ] = max(term0(:,2*tptest));
            T = term0(t0ix,(2*tptest)-1) - term0(t0in,(2*tptest)-1); % T = T(Emax) - T(root)

            % set phi IC
            optIC_four_file = [pwd '/data/spectrum/slingshots/optICfourfull_2048_LuFour_long_E0(1000)_Timept_3_lambda(' num2str(lambda) ').dat'];
            optIC_four = readmatrix(optIC_four_file);  
            solution = optIC_four( : , t0in );
            solution = solution';

            % Handle different resolutions
            if length(x) ~= length(x_ref)
                solution = real(ifft(solution));
                solution = interp1(x_ref, solution , x, 'spline');
                solution = fft(solution);
            end 

            %{
            optIC_phys_file = [pwd '/data/spectrum/slingshots/optICphys_' ...
                num2str(length(x_ref)*1) '_LuFour_long_E0(1000)_Timept_3_lambda(0.5).dat'];
            optIC_phys = readmatrix(optIC_phys_file);         
            optIC_four_file = [pwd '/data/spectrum/slingshots/optICfourfull_' ...
                num2str(length(x_ref)*1) '_LuFour_long_E0(' num2str(1000) ')_Timept_3_lambda(0.5).dat'];
            optIC_four = readmatrix(optIC_four_file);  

            % use end/4 [will not be < 1e-10 when growth starts] to find enstrophy root
            slingshot_index = find( abs(optIC_phys(end/4,2:end-1)) > 1e-10 , 1, 'first'); 
            solution = optIC_four( : , slingshot_index+1 );
            solution = solution';
            
            % Handle different resolutions
            if length(x) ~= length(x_ref)
                solution = real(ifft(solution));
                solution = interp1(x_ref, solution , x, 'spline');
                solution = fft(solution);
            end            

            % use branch data to find T = T(Emax) - T(root)
            maxbranch_file = [pwd '/data/enstrophy_solution/maxenstrophy_' num2str(length(x_ref)*1) '_LuFour_long_E0(' ...
                num2str(1000) ')_lambda(0.5).dat'];
            maxbranch = readmatrix(maxbranch_file); 
            finalbranch_file = [pwd '/data/enstrophy_solution/finalenstrophy_' num2str(length(x_ref)*1) '_LuFour_long_E0(' ...
                num2str(1000) ')_lambda(0.5).dat'];
            finalbranch = readmatrix(finalbranch_file); 
            TW = maxbranch(tptest,1);
            Tmax_index = maxbranch(tptest,4); % assuming stepscale = 1 (i.e. testcase has no '_tX')
            TW_pts = finalbranch(tptest,4); % number of timepoints in enstrophy evolution
            t_max = linspace(0,TW,TW_pts); % t for max index
            t_root = linspace(0,TW,size(optIC_phys,2)-4); % t for root index
            % subtract 3 columns from optIC_phys file
            Tmax = t_max(Tmax_index); 
            T0 = t_root((slingshot_index+1)-3); 
            T = Tmax - T0;
            %}
        case 'slingshot2'

            tptest = 10;
            % find indices of abs(min) and max pts in enstrophy evolution
            term0_file = [pwd '/data/time_evolution/terminal_2048_slingshotTlong_long3_E0(1000)_lambda(' num2str(lambda) ').dat'];  
            term0 = readmatrix(term0_file);
            [ ~ , t0in ] = min( abs( term0(3:end,2*tptest) ) );
            t0in = t0in+8;
            [ ~ , t0ix ] = max(term0(:,2*tptest));
            T = term0(t0ix,(2*tptest)-1) - term0(t0in,(2*tptest)-1); % T = T(Emax) - T(root)

            % set phi IC
            optIC_four_file = [pwd '/data/spectrum/slingshots/optICfourfull_2048_slingshotTlong_long3_E0(1000)_Timept_10_lambda(' num2str(lambda) ').dat'];
            optIC_four = readmatrix(optIC_four_file);  
            solution = optIC_four( : , t0in );
            solution = solution';

            % Handle different resolutions
            if length(x) ~= length(x_ref)
                solution = real(ifft(solution));
                solution = interp1(x_ref, solution , x, 'spline');
                solution = fft(solution);
            end 

            %{
            optIC_phys_file = [pwd '/data/spectrum/spectrum_E0_' num2str(1000) '/optICphys_' ...
                num2str(length(x_ref)*1) '_slingshotTlong_long3_E0(1000)_Timept_10_lambda(0.5).dat'];
            optIC_phys = readmatrix(optIC_phys_file);         
            optIC_four_file = [pwd '/data/spectrum/spectrum_E0_' num2str(1000) '/optICfourfull_' ...
                num2str(length(x_ref)*1) '_slingshotTlong_long3_E0(' num2str(1000) ')_Timept_10_lambda(0.5).dat'];
            optIC_four = readmatrix(optIC_four_file);  

            % use end/4 [will not be < 1e-10 when growth starts] to find enstrophy root
            slingshot_index = find( abs(optIC_phys(end/4,2:end-1)) > 1e-10 , 1, 'first'); 
            solution = optIC_four( : , slingshot_index+1 );
            solution = solution';
            
            % Handle different resolutions
            if length(x) ~= length(x_ref)
                solution = real(ifft(solution));
                solution = interp1(x_ref, solution , x, 'spline');
                solution = fft(solution);
            end 

            % use branch data to find T = T(Emax) - T(root)
            maxbranch_file = [pwd '/data/enstrophy_solution/maxenstrophy_' num2str(length(x_ref)*1) '_slingshotTlong_long3_E0(1000)_lambda(0.5).dat'];
            maxbranch = readmatrix(maxbranch_file); 
            finalbranch_file = [pwd '/data/enstrophy_solution/finalenstrophy_' num2str(length(x_ref)*1) '_slingshotTlong_long3_E0(' ...
                num2str(1000) ')_lambda(0.5).dat'];
            finalbranch = readmatrix(finalbranch_file); 
            TW = maxbranch(tptest,1);
            Tmax_index = maxbranch(tptest,4); % assuming stepscale = 1 (i.e. testcase has no '_tX')
            TW_pts = finalbranch(tptest,4); % number of timepoints in enstrophy evolution
            t_max = linspace(0,TW,TW_pts); % t for max index
            t_root = linspace(0,TW,size(optIC_phys,2)-4); % t for root index
            % subtract 3 columns from optIC_phys file
            Tmax = t_max(Tmax_index); 
            T0 = t_root((slingshot_index+1)-3); 
            T = Tmax - T0; 
            %}
        case 'slingshot3'

            tptest = 29;
            % find indices of abs(min) and max pts in enstrophy evolution
            term0_file = [pwd '/data/time_evolution/terminal_2048_slingshot2Tlong2_long10_E0(1000)_lambda(' num2str(lambda) ').dat'];  
            term0 = readmatrix(term0_file);
            [ ~ , t0in ] = min( abs( term0(3:end,2*tptest) ) );
            t0in = t0in+2;
            [ ~ , t0ix ] = max(term0(:,2*tptest));
            T = term0(t0ix,(2*tptest)-1) - term0(t0in,(2*tptest)-1); % T = T(Emax) - T(root)

            % set phi IC
            optIC_four_file = [pwd '/data/spectrum/slingshots/optICfourfull_2048_slingshot2Tlong2_long10_E0(1000)_Timept_29_lambda(' num2str(lambda) ').dat'];
            optIC_four = readmatrix(optIC_four_file);  
            solution = optIC_four( : , t0in );
            solution = solution';

            % Handle different resolutions
            if length(x) ~= length(x_ref)
                solution = real(ifft(solution));
                solution = interp1(x_ref, solution , x, 'spline');
                solution = fft(solution);
            end

            %{
            optIC_phys_file = [pwd '/data/spectrum/spectrum_E0_' num2str(1000) '/optICphys_' ...
                num2str(length(x_ref)*1) '_slingshot2Tlong2_long10_E0(1000)_Timept_29_lambda(0.5).dat'];
            optIC_phys = readmatrix(optIC_phys_file);         
            optIC_four_file = [pwd '/data/spectrum/spectrum_E0_' num2str(1000) '/optICfourfull_' ...
                num2str(length(x_ref)*1) '_slingshot2Tlong2_long10_E0(' num2str(1000) ')_Timept_29_lambda(0.5).dat'];
            optIC_four = readmatrix(optIC_four_file);  

            % use end/4 [will not be < 1e-10 when growth starts] to find enstrophy root
            slingshot_index = find( abs(optIC_phys(end/4,2:end-1)) > 1e-10 , 1, 'first'); 
            solution = optIC_four( : , slingshot_index+1 );
            solution = solution';
            
            % Handle different resolutions
            if length(x) ~= length(x_ref)
                solution = real(ifft(solution));
                solution = interp1(x_ref, solution , x, 'spline');
                solution = fft(solution);
            end 

            % use branch data to find T = T(Emax) - T(root)
            maxbranch_file = [pwd '/data/enstrophy_solution/maxenstrophy_' num2str(length(x_ref)*1) '_slingshot2Tlong2_long10_E0(' ...
                num2str(1000) ')_lambda(0.5).dat'];
            maxbranch = readmatrix(maxbranch_file); 
            finalbranch_file = [pwd '/data/enstrophy_solution/finalenstrophy_' num2str(length(x_ref)*1) '_slingshot2Tlong2_long10_E0(' ...
                num2str(1000) ')_lambda(0.5).dat'];
            finalbranch = readmatrix(finalbranch_file); 
            TW = maxbranch(tptest,1);
            Tmax_index = maxbranch(tptest,4); % assuming stepscale = 1 (i.e. testcase has no '_tX')
            TW_pts = finalbranch(tptest,4); % number of timepoints in enstrophy evolution
            t_max = linspace(0,TW,TW_pts); % t for max index
            t_root = linspace(0,TW,size(optIC_phys,2)-4); % t for root index
            % subtract 3 columns from optIC_phys file
            Tmax = t_max(Tmax_index); 
            T0 = t_root((slingshot_index+1)-3); 
            T = Tmax - T0;
            %}
        case 'slingshot4'

            tptest = 19;
            % find indices of abs(min) and max pts in enstrophy evolution
            term0_file = [pwd '/data/time_evolution/terminal_2048_slingshot3Tfull_long29_E0(1000)_lambda(' num2str(lambda) ').dat'];  
            term0 = readmatrix(term0_file);
            [ ~ , t0in ] = min( abs( term0(3:end,2*tptest) ) );
            t0in = t0in+14;
            [ ~ , t0ix ] = max(term0(:,2*tptest));
            T = term0(t0ix,(2*tptest)-1) - term0(t0in,(2*tptest)-1); % T = T(Emax) - T(root)

            % set phi IC
            optIC_four_file = [pwd '/data/spectrum/slingshots/optICfourfull_2048_slingshot3Tfull_long29_E0(1000)_Timept_19_lambda(' num2str(lambda) ').dat'];
            optIC_four = readmatrix(optIC_four_file);  
            solution = optIC_four( : , t0in );
            solution = solution';

            % Handle different resolutions
            if length(x) ~= length(x_ref)
                solution = real(ifft(solution));
                solution = interp1(x_ref, solution , x, 'spline');
                solution = fft(solution);
            end


            %{
            optIC_phys_file = [pwd '/data/spectrum/spectrum_E0_' num2str(1000) '/optICphys_' ...
                num2str(length(x_ref)*1) '_slingshot3Tfull_long29_E0(1000)_Timept_19_lambda(0.5).dat'];
            optIC_phys = readmatrix(optIC_phys_file);         
            optIC_four_file = [pwd '/data/spectrum/spectrum_E0_' num2str(1000) '/optICfourfull_' ...
                num2str(length(x_ref)*1) '_slingshot3Tfull_long29_E0(' num2str(1000) ')_Timept_19_lambda(0.5).dat'];
            optIC_four = readmatrix(optIC_four_file);  

            % use end/4 [will not be < 1e-10 when growth starts] to find enstrophy root
            slingshot_index = find( abs(optIC_phys(end/4,2:end-1)) > 1e-10 , 1, 'first'); 
            solution = optIC_four( : , slingshot_index+1 );
            solution = solution';
            
            % Handle different resolutions
            if length(x) ~= length(x_ref)
                solution = real(ifft(solution));
                solution = interp1(x_ref, solution , x, 'spline');
                solution = fft(solution);
            end 

            % use branch data to find T = T(Emax) - T(root)
            maxbranch_file = [pwd '/data/enstrophy_solution/maxenstrophy_' num2str(length(x_ref)*1) '_slingshot3Tfull_long29_E0(' ...
                num2str(1000) ')_lambda(0.5).dat'];
            maxbranch = readmatrix(maxbranch_file); 
            finalbranch_file = [pwd '/data/enstrophy_solution/finalenstrophy_' num2str(length(x_ref)*1) '_slingshot3Tfull_long29_E0(' ...
                num2str(1000) ')_lambda(0.5).dat'];
            finalbranch = readmatrix(finalbranch_file); 
            TW = maxbranch(tptest,1);
            Tmax_index = maxbranch(tptest,4); % assuming stepscale = 1 (i.e. testcase has no '_tX')
            TW_pts = finalbranch(tptest,4); % number of timepoints in enstrophy evolution
            t_max = linspace(0,TW,TW_pts); % t for max index
            t_root = linspace(0,TW,size(optIC_phys,2)-4); % t for root index
            % subtract 3 columns from optIC_phys file
            Tmax = t_max(Tmax_index); 
            T0 = t_root((slingshot_index+1)-3); 
            T = Tmax - T0;
            %}
        case 'slingshot5'


            tptest = 27;
            % find indices of abs(min) and max pts in enstrophy evolution
            term0_file = [pwd '/data/time_evolution/terminal_2048_slingshot4Tlong_full19_E0(1000)_lambda(' num2str(lambda) ').dat'];  
            term0 = readmatrix(term0_file);
            [ ~ , t0in ] = min( abs( term0(3:end,2*tptest) ) );
            t0in = t0in+2;
            [ ~ , t0ix ] = max(term0(:,2*tptest));
            T = term0(t0ix,(2*tptest)-1) - term0(t0in,(2*tptest)-1); % T = T(Emax) - T(root)

            % set phi IC
            optIC_four_file = [pwd '/data/spectrum/slingshots/optICfourfull_2048_slingshot4Tlong_full19_E0(1000)_Timept_27_lambda(' num2str(lambda) ').dat'];
            optIC_four = readmatrix(optIC_four_file);  
            solution = optIC_four( : , t0in );
            solution = solution';

            % Handle different resolutions
            if length(x) ~= length(x_ref)
                solution = real(ifft(solution));
                solution = interp1(x_ref, solution , x, 'spline');
                solution = fft(solution);
            end

            %{
            optIC_phys_file = [pwd '/data/spectrum/spectrum_E0_' num2str(1000) '/optICphys_' ...
                num2str(length(x_ref)*1) '_slingshot4Tlong_full19_E0(1000)_Timept_27_lambda(0.5).dat'];
            optIC_phys = readmatrix(optIC_phys_file);         
            optIC_four_file = [pwd '/data/spectrum/spectrum_E0_' num2str(1000) '/optICfourfull_' ...
                num2str(length(x_ref)*1) '_slingshot4Tlong_full19_E0(' num2str(1000) ')_Timept_27_lambda(0.5).dat'];
            optIC_four = readmatrix(optIC_four_file);  

            % use end/4 [will not be < 1e-10 when growth starts] to find enstrophy root
            slingshot_index = find( abs(optIC_phys(end/4,2:end-1)) > 1e-10 , 1, 'first'); 
            solution = optIC_four( : , slingshot_index+1 );
            solution = solution';
            
            % Handle different resolutions
            if length(x) ~= length(x_ref)
                solution = real(ifft(solution));
                solution = interp1(x_ref, solution , x, 'spline');
                solution = fft(solution);
            end 

            % use branch data to find T = T(Emax) - T(root)
            maxbranch_file = [pwd '/data/enstrophy_solution/maxenstrophy_' num2str(length(x_ref)*1) '_slingshot4Tlong_full19_E0(' ...
                num2str(1000) ')_lambda(0.5).dat'];
            maxbranch = readmatrix(maxbranch_file); 
            finalbranch_file = [pwd '/data/enstrophy_solution/finalenstrophy_' num2str(length(x_ref)*1) '_slingshot4Tlong_full19_E0(' ...
                num2str(1000) ')_lambda(0.5).dat'];
            finalbranch = readmatrix(finalbranch_file); 
            TW = maxbranch(tptest,1);
            Tmax_index = maxbranch(tptest,4); % assuming stepscale = 1 (i.e. testcase has no '_tX')
            TW_pts = finalbranch(tptest,4); % number of timepoints in enstrophy evolution
            t_max = linspace(0,TW,TW_pts); % t for max index
            t_root = linspace(0,TW,size(optIC_phys,2)-4); % t for root index
            % subtract 3 columns from optIC_phys file
            Tmax = t_max(Tmax_index); 
            T0 = t_root((slingshot_index+1)-3); 
            T = Tmax - T0;
            %}
        case 'slingshot6'

            tptest = 31;
            % find indices of abs(min) and max pts in enstrophy evolution
            term0_file = [pwd '/data/time_evolution/terminal_2048_slingshot4Tlong_full19_E0(1000)_lambda(' num2str(lambda) ').dat'];  
            term0 = readmatrix(term0_file);
            [ ~ , t0in ] = min( abs( term0(3:end,2*tptest) ) );
            t0in = t0in+2;
            [ ~ , t0ix ] = max(term0(:,2*tptest));
            T = term0(t0ix,(2*tptest)-1) - term0(t0in,(2*tptest)-1); % T = T(Emax) - T(root)

            % set phi IC
            optIC_four_file = [pwd '/data/spectrum/slingshots/optICfourfull_2048_slingshot4Tlong_full19_E0(1000)_Timept_31_lambda(' num2str(lambda) ').dat'];
            optIC_four = readmatrix(optIC_four_file);  
            solution = optIC_four( : , t0in );
            solution = solution';

            % Handle different resolutions
            if length(x) ~= length(x_ref)
                solution = real(ifft(solution));
                solution = interp1(x_ref, solution , x, 'spline');
                solution = fft(solution);
            end


            %{
            optIC_phys_file = [pwd '/data/spectrum/spectrum_E0_' num2str(1000) '/optICphys_' ...
                num2str(length(x_ref)*1) '_slingshot4Tlong_full19_E0(1000)_Timept_31_lambda(0.5).dat'];
            optIC_phys = readmatrix(optIC_phys_file);         
            optIC_four_file = [pwd '/data/spectrum/spectrum_E0_' num2str(1000) '/optICfourfull_' ...
                num2str(length(x_ref)*1) '_slingshot4Tlong_full19_E0(' num2str(1000) ')_Timept_31_lambda(0.5).dat'];
            optIC_four = readmatrix(optIC_four_file);  

            % use end/4 [will not be < 1e-10 when growth starts] to find enstrophy root
            slingshot_index = find( abs(optIC_phys(end/4,2:end-1)) > 1e-10 , 1, 'first'); 
            solution = optIC_four( : , slingshot_index+1 );
            solution = solution';
           
            % Handle different resolutions
            if length(x) ~= length(x_ref)
                solution = real(ifft(solution));
                solution = interp1(x_ref, solution , x, 'spline');
                solution = fft(solution);
            end 

            % use branch data to find T = T(Emax) - T(root)
            maxbranch_file = [pwd '/data/enstrophy_solution/maxenstrophy_' num2str(length(x_ref)*1) '_slingshot4Tlong_full19_E0(' ...
                num2str(1000) ')_lambda(0.5).dat'];
            maxbranch = readmatrix(maxbranch_file); 
            finalbranch_file = [pwd '/data/enstrophy_solution/finalenstrophy_' num2str(length(x_ref)*1) '_slingshot4Tlong_full19_E0(' ...
                num2str(1000) ')_lambda(0.5).dat'];
            finalbranch = readmatrix(finalbranch_file); 
            TW = maxbranch(tptest,1);
            Tmax_index = maxbranch(tptest,4); % assuming stepscale = 1 (i.e. testcase has no '_tX')
            TW_pts = finalbranch(tptest,4); % number of timepoints in enstrophy evolution
            t_max = linspace(0,TW,TW_pts); % t for max index
            t_root = linspace(0,TW,size(optIC_phys,2)-4); % t for root index
            % subtract 3 columns from optIC_phys file
            Tmax = t_max(Tmax_index); 
            T0 = t_root((slingshot_index+1)-3); 
            T = Tmax - T0;
            %}
    end
return