function main_1DBEnstrophy_2024_v3

%{
[v3: clean up code and documentation]
Update: Jovan Zigic, 2024/05
Contributors: Jovan Zigic, Diego Ayala, Bartosz Protas
Description: 
    Computing Enstrophy Growth (subject to initial enstrophy condition) via 
    Riemannian Gradient Optimization of Heuristic Choice of Initial Data for 
    1D viscous Burgers (in Fourier domain).
Papers:
    Maximum Enstrophy Growth in Burgers Equation, D. Ayala, 2010 (MS Thesis)
    On maximum enstrophy growth in a hydrodynamic system, D. Ayala and B. Protas, 2011 (Physica D)
Design (by level):
1. Choose initial data s, shift time window, and set smoothing parameter.
    Default: s = 0 <-> initial data from (Doering, Lu 2009).
    Default: shiftsign = 0 <-> no time window shift.
    Default: lamtest = 3 <-> lambda = 1.
2. Choose whether to shift time window T or test slingshot points.
    Default: shiftT = 0 || tptest = 1 <-> not shifting or slingshot testing.
    Tuning parameters for slingshot testing: 
        {lamtest, tptest, ens, timept,
        ss_shift, TimeWindow # of points, testcase, CONT, shotname}
3. Choose initial enstrophy values for computation (or shift slingshot initial data).
    Default: E0 = 1000 <-> compute for initial enstrophy of 1000.
    Default: ss_shift = 2 <-> no shift in slingshot initial data.
4. Choose branch of time windows for E0 level.
    Default: T = [2/(50*sqrt(E0)), 2/(1*sqrt(E0))] <-> maximum enstrophy for E0 expected at T ~ 1/sqrt(E0).
5. Choose continuation or no continuation of initial data for branch computation.
    Default: p = 1:2 <-> p=1 uses previous timept solution, p=2 uses (Doering, Lu 2009) initial data.
6. Search for optimal solution using gradient optimization routine
%}

    s = 0; % slingshot number (0 if not using slingshot initial data)
    shiftsign = 0; % choose whether to shift time window (0: no shift, 1: positive shift, 2: negative shift)
    % lamtest = 3; % choice of smoothing parameter

    for lamtest = 3:3 % choice of smoothing parameter

        % sensitivity testing of time window perturbation % 
        if shiftsign < 2
            pm = 'p'; % positive shift
        else
            pm = 'm';  % negative shift
        end

        % shiftT = 0; % time window shift parameter
        tptest = 1; % slingshot search parameter

        for shiftT = 0:0 % "shiftT" <-> sensitivity testing || "tptest" <-> slingshot searching
            clearvars -except shiftT tptest s shiftsign pm lamtest; 
            close all;

            % initialize diagnostic switches %
            Nscale = 1; % physical domain resolution scaling
            stepscale = 1; % 1/stepscale = time-step scaling for oscillations in enstrophy values
            VERBOSE = 1; % set level of data output

            % declare and initialize parameters %
            CONTSET = 0; % ( 0: Lu, 2: other data file, 3: slingshot testing )
            ensstart = 1;
            ensend = 1;
            timestart = 1; % 
            timeend = 3*19;% 

            N = 2048*Nscale; % number of grid points
            MAXITER = 1000; % this should never be reached
            J_step_tol = 10^(-6); % 10^(-6) is sufficient
            x = linspace(0,1-1/N,N); % physical space domain
            x_2048 = linspace(0,1-1/2048,2048); % required for LuLu optIC interpolation to higher N
            K0 = [0:N/2-1 0 -N/2+1:-1]; % fourier space domain
            K1 = [0:N/2 -N/2+1:-1]; % fourier space domain for derivative 1
            nu = 0.001; % viscosity coefficient
            lam = [ 2, 1.5, 1, 0.5, 0.25, 0.1 ];
            lambda = lam(lamtest); % smoothing parameter
            % E0 = % initial enstrophy
            % T = % length of time window
            % uField = % solution to unknown variable
            % tvector = % solution timepoints

            EnstPoints = 31; % number of initial enstrophy points
            Enstrophy = logspace(3,6,EnstPoints);
                            
            % turn on/off slingshot initial data choice %
            if CONTSET < 3
                ss_shift = 0;
            else
                ss_shift = 2; 
            end
            % enspt = ensstart; % choose enstrophy point to test

            for enspt = ensstart:ensend % "enspt" <-> enstrophy point || "ss_shift" <-> shift slingshot initial data

                % choose time windows to test at each E0 %
                E0 = Enstrophy(enspt);
                prefactor = 2;
                TimePoints = 50; % number of timepoints in each initial enstrophy branch
                T_ens_UB = prefactor*(1/sqrt(E0)); % use for max enstrophy T 
                T_ens_LB = T_ens_UB/TimePoints; % Adjust for decreasing T_max  
                TimeWindow = linspace(T_ens_LB,T_ens_UB,3*TimePoints); 
        
                % choose data output at each E0 %
                switch VERBOSE
                    case 1 
                        mkdir([pwd  '/data/diagnostics/diagnostics_E0_' num2str(E0) '' ]);
                        mkdir([pwd  '/data/enstrophy_solution' ]);
                        mkdir([pwd  '/data/time_evolution' ]);
                        mkdir([pwd  '/data/spectrum' ]);
                        mkdir([pwd  '/data/spectrum/spectrum_E0_' num2str(E0) '' ]);
                        mkdir([pwd  '/data/spectrum/temp' ]);
                        mkdir([pwd  '/data/kappa/kappa_E0_' num2str(E0) '' ]);
                        mkdir([pwd  '/data/runtime' ]);
                    case 3
                        mkdir([pwd  '/data/kappa' ]);
                end
        
                for timept = timestart:timeend 

                    ParamPoints = 2; 
                    enstest = NaN(1,ParamPoints);
                    phitest = NaN(ParamPoints,N);

                    %%% name testcase for output data %%%
                    if s == 0
                        testcase_max = [ '_' num2str(N) '_LuMax_t' num2str(stepscale) '_new' ];
                        testcase1 = [ '_' num2str(N) '_LuCont1_t' num2str(stepscale) '_new' ];
                        testcase0 = [ '_' num2str(N) '_LuCont0_t' num2str(stepscale) '_new' ];
                        %%% other initial data if CONTSET == 2 %%%
                        testcase1_phi = [ '_' num2str(N) '_LuCont1_t' num2str(stepscale) '_short120' ];
                        testcase0_phi = [ '_' num2str(N) '_LuCont0_t' num2str(stepscale) '_short120' ];
                        %{
                        %%% slingshot testing %%%
                        testcase_max = [ '_' num2str(N) '_L1Max_t' num2str(stepscale) '_tp(' num2str(tptest) ')' ];
                        testcase1 = [ '_' num2str(N) '_L1Cont1_t' num2str(stepscale) '_tp(' num2str(tptest) ')' ]; 
                        testcase0 = [ '_' num2str(N) '_L1Cont0_t' num2str(stepscale) '_tp(' num2str(tptest) ')' ];
                        %}
                        if ss_shift ~= 0
                            testcase_max = [ '_' num2str(N) '_L1Max_t' num2str(stepscale) '_tp(' num2str(tptest) ')_shift' num2str(ss_shift) '' ];
                            testcase1 = [ '_' num2str(N) '_L1Cont1_t' num2str(stepscale) '_tp(' num2str(tptest) ')_shift' num2str(ss_shift) '' ]; 
                            testcase0 = [ '_' num2str(N) '_L1Cont0_t' num2str(stepscale) '_tp(' num2str(tptest) ')_shift' num2str(ss_shift) '' ];
                        end
                        if shiftT ~= 0
                            testcase_max = [ '_' num2str(N) '_LuMax_t' num2str(stepscale) '_shift' pm '' num2str(shiftT) '' ];
                            testcase1 = [ '_' num2str(N) '_LuMax_t' num2str(stepscale) '_shift' pm '' num2str(shiftT) '' ];
                            testcase0 = [ '_' num2str(N) '_LuMax_t' num2str(stepscale) '_shift' pm '' num2str(shiftT) '' ];
                        end
                    else
                        testcase_max = [ '_' num2str(N) '_slingshot' num2str(s) 'Max_t' num2str(stepscale) '_fine50' ];
                        testcase1 = [ '_' num2str(N) '_slingshot' num2str(s) 'Cont1_t' num2str(stepscale) '_fine50' ];
                        testcase0 = [ '_' num2str(N) '_slingshot' num2str(s) 'Cont0_t' num2str(stepscale) '_fine50' ];
                        testcase1_phi = [ '_' num2str(N) '_slingshot' num2str(s) 'Cont1_t' num2str(stepscale) '_short120' ];
                        testcase0_phi = [ '_' num2str(N) '_slingshot' num2str(s) 'Cont0_t' num2str(stepscale) '_short120' ];
                        if shiftT ~= 0
                            testcase_max = [ '_' num2str(N) '_slingshot' num2str(s) 'Max_t' num2str(stepscale) ...
                                '_shift' pm '' num2str(shiftT) '' ];
                            testcase1 = [ '_' num2str(N) '_slingshot' num2str(s) 'Cont1_t' num2str(stepscale) ...
                                '_shift' pm '' num2str(shiftT) '' ];
                            testcase0 = [ '_' num2str(N) '_slingshot' num2str(s) 'Cont0_t' num2str(stepscale) ...
                                '_shift' pm '' num2str(shiftT) '' ];
                        end
                    end

                    for p = 1:ParamPoints
                
                        if p == 1
                            CONT = 1;
                            testcase = testcase1;
                            testcase_phi = testcase1_phi;
                        elseif p == 2
                            CONT = CONTSET;
                            testcase = testcase0;
                            testcase_phi = testcase0_phi;
                        end
    
                        tic % start timer
                        % initialize diagnostics
                        J = NaN(MAXITER);
                        K = NaN(MAXITER);
                        step = NaN(MAXITER);
                        momentum = NaN(MAXITER);
                        
                        % Initial condition of physical problem %
                        if ( timept == timestart || CONT ~= 1 ) 
                            if ( CONT == 0 || CONT == 1 )
                                [ phi , ~ ] = initialguess('exact', x, E0, 0, x_2048, timept,tptest,testcase,lambda,ss_shift);
                            elseif CONT == 2
                                [ phi , ~ ] = initialguess('optIC', x, E0, 0, x_2048, timept,tptest,testcase_phi,lambda,ss_shift);
                            elseif CONT == 3 
                                %%% choose slingshot setting %%%
                                % shotname = ['slingshot' num2str(s) ''];
                                shotname = 'ss_shift';
                                % shotname = 'ss_root';
                                % shotname = 'ss_start';
                                % shotname = 'slingshot';
                                [ phi , T_interval ] = initialguess(shotname, x, E0, 0, x_2048, timept,tptest,testcase,lambda,ss_shift);
                                % adjust time interval only when testing slingshot points
                                if timept == timestart
                                    TimeWindow = linspace( T_interval/2 , 2*T_interval , timeend ); 
                                end
                            end
    
                            phi = adjust_optIC(phi,E0,K0,N); % initial data 
                            phi_x = 2*pi*1i*K0.*phi; % derivative of initial data
                            E = 0.5*sum(abs(phi_x).^2)/N^2; % enstrophy E0 for retraction
                        end
    
                        if shiftsign < 2
                            T = TimeWindow(timept)*(1000+shiftT)/1000;
                        else
                            T = TimeWindow(timept)*(1000-shiftT)/1000;
                        end

                        ITER = 1; % start counting iterations
                        
                        fprintf('\n\n  Enstrophy point = %d,  Time point = %d \n',enspt,timept);
                        fprintf('  Iterations: ')
    
                        % Do at least one iteration
                        [tvector,uField] = BurgersDS_Fourier(phi,K1,K0,T,nu,N,stepscale); % solve PDE
                        Ntime = length(tvector); % number of time-steps in PDE solution 
                        uu = uField(Ntime,:); % physical solution at final time-step
                        J(ITER) = eval_J(uu,phi,K0,N); % evaluate objective functional (E(T) - E0)
                        K(ITER) = norm(uu); % kinetic energy at E(T)
                        step(ITER) = 0;
                        momentum(ITER) = 0;
                        
                        u_adj = BurgersAS_Fourier(uu,K1,K0,T,nu,N,uField,tvector,stepscale); % solve adjoint PDE
                        gradJ = eval_grad_J(u_adj,phi,K1,lambda); % gradient of objective functional
                        gradJ_x = (2*pi*1i*K0).*gradJ; % derivative of gradient of objective functional
                        
                        %%% Retracted Conjugate Gradient algorithm %%%
                        phi_x = 2*pi*1i*K0.*phi; % derivative of physical data
                        tau = -2*real( sum(phi.*conj(gradJ)) + lambda^2*sum(phi_x.*conj(gradJ_x)) ) ... 
                            /( sum(abs(gradJ).^2) + lambda^2*sum(abs(gradJ_x).^2) ) ; % initial guess for optimal step length
                        direction = gradJ; % direction towards optimal solution
                        tau = optimize_tau(phi,direction,abs(0.2*tau),E,K1,K0,T,nu,N,stepscale); % compute initial optimal step size 
                        phi = phi + tau*direction; % iterative solution update
                        phi = adjust_optIC(phi,E,K0,N); % retraction to constraint manifold
                        
                        [tvector,uField] = BurgersDS_Fourier(phi,K1,K0,T,nu,N,stepscale);
                        Ntime = length(tvector);
                        uu = uField(Ntime,:);
                        J_new = eval_J(uu, phi,K0,N);
                        
                        J_step = (J_new - J(ITER))/J(ITER); % relative change in objective functional value
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
                                % test definition of Gateaux differential
                                [epsilon,kappa] = kappaTestFourier(phi,gradJ,gradJ_x,J(ITER),x,K1,K0,T,nu,N,lambda,stepscale);
                                kappa_save(kappa, epsilon, timept,E0,lambda,ITER,testcase);
                        end
                    
                        while abs(J_step) > J_step_tol && ITER <= MAXITER
                            
                            alpha = sum(abs(gradJ).^2)/N^2 + lambda^2*sum(abs(gradJ_x).^2)/N^2; % normalizer for beta
                            u_adj = BurgersAS_Fourier(uu,K1,K0,T,nu,N,uField,tvector,stepscale);
                            gradJ_new = eval_grad_J(u_adj,phi,K1,lambda);
                            gradJ_x_new = (2*pi*1i*K0).*gradJ_new;
                            beta = sum(gradJ_new.*conj(gradJ_new - gradJ))/N^2 ...
                                + lambda^2*sum(gradJ_x_new.*conj(gradJ_x_new - gradJ_x))/N^2; % Polak-Ribiere conjugate direction
                            gradJ = gradJ_new;
                            gradJ_x = gradJ_x_new;
                            
                            switch VERBOSE
                                case 2 % currently out of order
                                    if ITER < 3
                                        [epsilon,kappa] = kappaTestFourier(phi,gradJ,gradJ_x,J(ITER),x,K1,K0,T,nu,N,lambda,stepscale);
                                        kappa_save(kappa, epsilon, timept,E0,lambda,ITER,testcase);
                                    end
                            end
                            
                            % reset beta after 20 iterations %
                            if mod(ITER,20)==0
                                beta=0;
                            else
                                beta = max(real(beta/alpha), 0);
                            end
                            
                            %%% Retracted Conjugate Gradient algorithm %%%
                            direction = gradJ - beta*direction; % direction towards optimal solution
                            tau = optimize_tau(phi,direction,tau,E,K1,K0,T,nu,N,stepscale); % compute optimal step size 
                            phi = phi + tau*direction; % iterative solution update
                            phi = adjust_optIC(phi,E,K0,N); % retraction to constraint manifold
                            
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
        
                            disp(['Completed iteration ' num2str(ITER-1)]);
                            disp(['Runtime: ' num2str(toc)]);
                        end
        
                        runtime = toc; % stop timer
                        
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
                            case 1 % save data %
                                diagnostics_save(J,K,step,momentum,timept,lambda,E0,ITER,testcase);
                                f_ens = time_evolution_and_optIC_save(Ntime,phi,timept,E0,lambda,T,N,uField,testcase,K0,stepscale);
                                branch_save(timept,f_ens,E0,lambda,T,testcase); 
                        end
            
                        phi_save(phi, timept, E0, timestart, testcase); % save optimal initial data 
                        enstest(1,p) = max(f_ens); % save max enstrophy value
                        phitest(p,:) = phi; % save optimal initial data for continuation
                        runtime_save(runtime, testcase, E0, timept, lambda, timestart); % save runtime data for each routine

                    end %parampt

                    % choose higher enstrophy result %
                    if enstest(1,1) > enstest(1,2) || isnan(enstest(1,2))
                        pc = 1;
                        testcase_old = testcase1;
                        testcase_other = testcase0;
                    else
                        pc = 2;
                        testcase_old = testcase0;
                        testcase_other = testcase1;
                    end

                    % save data with higher enstrophy (continuation or no continuation) %
                    phi = phitest(pc,:);
                    max_save(timept,lambda,E0,testcase_max,testcase_old,testcase_other,timeend); 

                    if timept == timeend
                        growth_save(lambda,E0,testcase_max); % save max enstrophy across branch of time windows   
                    end

                end %timept
            end %enspt
        end 
    end
return