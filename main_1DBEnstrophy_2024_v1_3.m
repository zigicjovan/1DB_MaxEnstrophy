function main_1DBEnstrophy_2024_v1_3

% v1_1: slingshot tests within specific time window
% v1_2: fixed physical solution output (time_evol*)
% v1_3: create smoothing test (shift T), fix overwrite of evolution file

    for ss = 0:0 % slingshot testing (0: original, 1: slingshot 1, etc. )
        if ss < 7
            s = ss;
            pm = 'p';
        else
            s = ss - 6;
            pm = 'm';
        end
        for shiftT = 0:0 % sensitivity testing
            clearvars -except shiftT s ss pm; 
            close all;

            % initialize diagnostic switches %
            Nscale = 1;
            stepscale = 1; % 1/stepscale = time-step scaling for oscillations in enstrophy values
            tptest = 3;
            VERBOSE = 1; 

            % declare and initialize parameters %
            CONT = 3; % no continuation recommended (0: Lu, 1: cont, 2: 2048 opt IC, 3: slingshot)

            if s == 2
                CONT = 0;
            elseif s == 0
                CONT = 0;
                stepscale = 1;
            end
            %}
            ensstart = 1;
            ensend = 1;
            timestart = 13; % max occurs ~ timepoint 11
            timeend = 13; % 31 for full/long, 15 for level test
            %
            if s == 1
                timestart = 8;
                timeend = 12;
            elseif s == 2
                timestart = 7 ;
                timeend = 14;
            elseif s == 3
                timestart = 8;
                timeend = 10;
            elseif s == 4
                timestart = 6;
                timeend = 16;
            elseif s == 5
                timestart = 6;
                timeend = 10;
            elseif s == 6
                timestart = 7;
                timeend = 15;
            end
            %}
            N = 2048*Nscale; % number of grid points
            MAXITER = 6; % 6 iterations is sufficient
            J_step_tol = 10^(-6); % 10^(-6) is sufficient
            x = linspace(0,1-1/N,N); % physical space domain
            x_2048 = linspace(0,1-1/2048,2048); % required for LuLu optIC interpolation to higher N
            K0 = [0:N/2-1 0 -N/2+1:-1]; % fourier space domain
            K1 = [0:N/2 -N/2+1:-1]; % fourier space domain for derivative 1
            nu = 0.001; % viscosity coefficient
            lambda = 0.5; % length-scale parameter
            % E0 = % initial enstrophy
            % T = % length of time window
            % uField = % solution to unknown variable
            % tvector = % solution timepoints

            % testcase = [ '_' num2str(N) '_slingshot5Tfull_long' num2str(tptest) '' ]; % " _[name] " % adjust prior to run
            % [slingshot] tpts = 31;  Tfull: cancel specific TW, Tlong: also set UB = 1
            if shiftT ~= 0
                testcase = [ '_' num2str(N) '_slingshot' num2str(s) ...
                    '_shift' pm '0' num2str(shiftT) '' ];
            end

            if Nscale > 1 || stepscale > 1
                testcase = [ '_' num2str(N) '_slingshot' num2str(s) '_t' num2str(stepscale) '' ];
                if shiftT ~= 0
                    testcase = [ '_' num2str(N) '_slingshot' num2str(s) '_t' num2str(stepscale) ...
                        '_shift' pm '' num2str(shiftT) '' ];
                end
            end
            if s == 0
                testcase = [ '_' num2str(N) '_Lu_t' num2str(stepscale) '' ];
            end
            
            ParamPoints = 1;
              
            for p = 1:ParamPoints
        
                % adjust what E0 is being tested
                if p == 1
                    EnstPoints = 31;
                    Enstrophy = logspace(3,6,EnstPoints);
                    ens_start = ensstart;
                    EnstPoints = ensend;
                elseif p == 2
                    EnstPoints = 31;
                    Enstrophy = logspace(3,6,EnstPoints);
                    ens_start = 16;
                    EnstPoints = 16;
                end
                
                for enspt = ens_start:EnstPoints
                    
                    E0 = Enstrophy(enspt);
                    prefactor = 2;
                    % number of timepoints in each initial enstrophy branch
                    TimePoints = 31;
                    T_ens_UB = prefactor*(1/sqrt(E0)); % use for max enstrophy T 
                    % T_ens_UB = 1; % use for longer T
                    T_ens_LB = T_ens_UB/TimePoints; % Adjust for decreasing T_max  
                    TimeWindow = linspace(T_ens_LB,T_ens_UB,TimePoints); 
                    % adjust what timepoints are being tested
                    time_start = timestart; 
                    TimePoints = timeend; 
            
                    switch VERBOSE
                        case 1 
                            mkdir([pwd  '/data/diagnostics/diagnostics_E0_' num2str(E0) '' ]);
                            mkdir([pwd  '/data/enstrophy_solution' ]);
                            mkdir([pwd  '/data/time_evolution' ]);
                            mkdir([pwd  '/data/spectrum' ]);
                            mkdir([pwd  '/data/spectrum/spectrum_E0_' num2str(E0) '' ]);
                            mkdir([pwd  '/data/kappa/kappa_E0_' num2str(E0) '' ]);
                            mkdir([pwd  '/data/runtime' ]);
                        case 3
                            mkdir([pwd  '/data/kappa' ]);
                    end
            
                    for timept = time_start:TimePoints % 
    
                        % manual continutation
                        %{
                        if enspt == 1
                            if timept > 25
                                CONT = 1;
                            elseif timept < 7
                                CONT = 1;
                            else
                                CONT = 3;
                            end
                        end
                        %}
    
                        tic
                        % initialize diagnostics
                        J = NaN(MAXITER);
                        K = NaN(MAXITER);
                        step = NaN(MAXITER);
                        momentum = NaN(MAXITER);
                        
                        % Initial condition of physical problem %
                        if ( timept == time_start || CONT ~= 1 )    
                            if CONT == 0 
                                [ phi , ~ ] = initialguess('exact', x, E0, 0, x_2048, timept,tptest);
                            elseif CONT == 2
                                [ phi , ~ ] = initialguess('optIC', x, E0, 0, x_2048, timept,tptest);
                            elseif CONT == 3
                                shotname = ['slingshot' num2str(s) ''];
                                [ phi , T_interval ] = initialguess(shotname, x, E0, 0, x_2048, timept,tptest);
                                % Adjust time interval
                                if timept == time_start
                                    % TimeWindow = linspace( T_interval/2 , 2*T_interval , timeend ); %specific TW
                                end
                            end

                            phi = adjust_optIC(phi,E0,K0,N);
                            phi_x = 2*pi*1i*K0.*phi;
                            E = 0.5*sum(abs(phi_x).^2)/N^2;
                        end
    
                        T = TimeWindow(timept);
                        if ss < 7
                            s = ss;
                            T = TimeWindow(timept)*(1000+shiftT)/1000;
                        else
                            s = ss - 6;
                            T = TimeWindow(timept)*(1000-shiftT)/1000;
                        end
                        ITER = 1;
                        
                        fprintf('\n\n  Enstrophy point = %d,  Time point = %d \n',enspt,timept);
                        fprintf('  Iterations: ')

                        % Do at least one iteration
                        [tvector,uField] = BurgersDS_Fourier(phi,K1,K0,T,nu,N,stepscale);
                        Ntime = length(tvector);
                        uu = uField(Ntime,:);
                        J(ITER) = eval_J(uu,phi,K0,N);
                        K(ITER) = norm(uu);
                        step(ITER) = 0;
                        momentum(ITER) = 0;
                        
                        u_adj = BurgersAS_Fourier(uu,K1,K0,T,nu,N,uField,tvector,stepscale);
                        gradJ = eval_grad_J(u_adj,phi,K1,lambda);
                        gradJ_x = (2*pi*1i*K0).*gradJ;
                        
                        phi_x = 2*pi*1i*K0.*phi;
                        tau = -2*real( sum(phi.*conj(gradJ)) + lambda^2*sum(phi_x.*conj(gradJ_x)) ) ... 
                            /( sum(abs(gradJ).^2) + lambda^2*sum(abs(gradJ_x).^2) ) ;
                        direction = gradJ;
                        tau = optimize_tau(phi,direction,abs(0.2*tau),E,K1,K0,T,nu,N,stepscale);
            
                        phi = phi + tau*direction;
                        phi = adjust_optIC(phi,E,K0,N);
                        
                        [tvector,uField] = BurgersDS_Fourier(phi,K1,K0,T,nu,N,stepscale);
                        Ntime = length(tvector);
                        uu = uField(Ntime,:);
                        J_new = eval_J(uu, phi,K0,N);
                        
                        J_step = (J_new - J(ITER))/J(ITER);
                        ITER = ITER + 1;
                        J(ITER) = J_new; 
                        K(ITER) = norm(uu);
                        step(ITER) = tau;
                        momentum(ITER) = 0;
        
                        switch VERBOSE
                            case 3
                                u_adj = BurgersAS_Fourier(uu,K1,K0,T,nu,N,uField,tvector,stepscale);
                                gradJ_new = eval_grad_J(u_adj,phi,K1,lambda);
                                gradJ_x_new = (2*pi*1i*K0).*gradJ_new;
                                gradJ = gradJ_new;
                                gradJ_x = gradJ_x_new;
                                [epsilon,kappa] = kappaTestFourier(phi,gradJ,gradJ_x,J(ITER),x,K1,K0,T,nu,N,lambda,stepscale);
                                kappa_save(kappa, epsilon, timept,E0,lambda,ITER,testcase);
                        end
                    
                        while abs(J_step) > J_step_tol && ITER <= MAXITER
                            
                            alpha = sum(abs(gradJ).^2)/N^2 + lambda^2*sum(abs(gradJ_x).^2)/N^2;
                            u_adj = BurgersAS_Fourier(uu,K1,K0,T,nu,N,uField,tvector,stepscale);
                            gradJ_new = eval_grad_J(u_adj,phi,K1,lambda);
                            gradJ_x_new = (2*pi*1i*K0).*gradJ_new;
                            beta = sum(gradJ_new.*conj(gradJ_new - gradJ))/N^2 ...
                                + lambda^2*sum(gradJ_x_new.*conj(gradJ_x_new - gradJ_x))/N^2;
                            gradJ = gradJ_new;
                            gradJ_x = gradJ_x_new;
                            
                            switch VERBOSE
                                case 1
                                    if ITER < 3
                                        [epsilon,kappa] = kappaTestFourier(phi,gradJ,gradJ_x,J(ITER),x,K1,K0,T,nu,N,lambda,stepscale);
                                        kappa_save(kappa, epsilon, timept,E0,lambda,ITER,testcase);
                                    end
                            end
                            
                            if mod(ITER,20)==0
                                beta=0;
                            else
                                beta = max(real(beta/alpha), 0);
                            end
                            
                            direction = gradJ - beta*direction;
                            tau = optimize_tau(phi,direction,tau,E,K1,K0,T,nu,N,stepscale);
                            
                            phi = phi + tau*direction;
                            phi = adjust_optIC(phi,E,K0,N);
                            
                            [tvector,uField] = BurgersDS_Fourier(phi,K1,K0,T,nu,N,stepscale);
                            Ntime = length(tvector);
                            uu = uField(Ntime,:);
                            J_new = eval_J(uu, phi,K0,N);
                            
                            J_step = (J_new - J(ITER))/J(ITER);
                            ITER = ITER+1;
                            J(ITER) = J_new;
                            K(ITER) = norm(uu);
                            step(ITER) = tau;
                            momentum(ITER) = beta;
            
                            if mod(ITER,50)==0
                                fprintf(' %d, ', ITER);                    
                            end
        
                            disp(['Completed iteration ' num2str(ITER-1)]);
                            disp(['Runtime: ' num2str(toc)]);
                        end
        
                        runtime = toc;
                        
                        disp('');
                        disp(['Time window: ' num2str(T)]);
                        disp(['Initial Enstrophy: ' num2str(E0)]);
                        disp(['Number of iterations: ' num2str(ITER-1)]);
                        disp(['Runtime: ' num2str(runtime)]);
                        
                        if ITER > MAXITER
                            disp('Optimization reached maximum number of iterations.');
                            disp(['Final Relative Step: ' num2str(J_step)]);
                        end
            
                        switch VERBOSE
                            case 1 
                                diagnostics_save(J,K,step,momentum,timept,lambda,E0,ITER,testcase);
                                f_ens = time_evolution_and_optIC_save(Ntime,phi,timept,E0,lambda,T,N,uField,testcase,K0,stepscale);
                                branch_save(timept,f_ens,E0,lambda,T,testcase); 
                        end
            
                        phi_save(phi, timept, E0, time_start, testcase);
                        runtime_save(runtime, testcase, E0, timept, lambda, time_start);
                    end %timepts
                    growth_save(lambda,E0,testcase);  
                end %enspts
            end %parampts
        end %ttest
    end
return