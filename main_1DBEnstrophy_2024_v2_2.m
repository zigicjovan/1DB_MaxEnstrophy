function main_1DBEnstrophy_2024_v2_2

%{
Updated by Jovan Zigic, 2024
Description (by loop level):
1. Choose initial data, and add +6 if shift in next level will be negative.
    Default: ss = 0 <-> initial data from (Doering, Lu 2009).
2. Choose whether to shift time window T or test slingshot points.
    Default: tptest = 1 || shiftT = 0 <-> not shifting or slingshot testing.
3. Choose initial enstrophy values for computation.
    Default: E0 = 1000 <-> compute for initial enstrophy of 1000.
4. Choose branch of time windows for computation for E0 level.
    Default: T = [2/(31*sqrt(E0)), 62/(31*sqrt(E0))] <-> maximum enstrophy for E0 expected at T ~ 1/sqrt(E0).
5. Choose continuation or no continuation of initial data for branch computation.
    Default: p = 1:2 <-> p=1 uses previous timept solution, p=2 uses (Doering, Lu 2009) initial data.
6. Search for optimal solution using method from (Ayala, Protas 2011).
%}

%{
v1_1: slingshot tests within specific time window
v1_2: fixed physical solution output (time_evol*)
v1_3: create smoothing test (shift T), fix overwrite of evolution file
v1_4: compute derivative (time_evol*), save max (cont/no cont)
v2_0: slingshot search version 2, temp spectrum files (max*, time*)
v2_1: evolution-spectrum match, delete Cont1 and Cont0 files (max_save)
v2_2: identify redundancies (slingshots, tptest), keep sensitivity testing
%}

%{
Current parameters:
lamtest
tptest
ens
timept
ss_shift
TimeWindow points
testcase
IC Cont
shotname
%}

    for lamtest = 4:5 % slingshot initial data (0: original (LuLu), 1: slingshot 1, etc.)
        ss = 0;
        if ss < 7
            s = ss; % slingshot number
            pm = 'p'; % positive shift
        else
            s = ss - 6; % slingshot number
            pm = 'm';  % negative shift
        end

        % tptest = 3; % slingshot search parameter
        shiftT = 0; % time window shift parameter

        for tptest = 19:19 % "shiftT" <-> sensitivity testing || "tptest" <-> slingshot searching
            clearvars -except ss_shift shiftT tptest s ss pm lamtest; 
            close all;

            % initialize diagnostic switches %
            Nscale = 1;
            stepscale = 1; % 1/stepscale = time-step scaling for oscillations in enstrophy values
            VERBOSE = 1; 

            % declare and initialize parameters %
            % [slingshot] tpts = 31;  Tfull: cancel specific TW, Tlong: also set UB = 1
            % CONT = 3; % ( 0: Lu, 1: cont, 2: 2048 opt IC, 3: slingshot )
            ensstart = 1;
            ensend = 1;
            timestart = 1; % 
            timeend = 40;% 
            %{
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
            MAXITER = 1000; % this should never be reached
            J_step_tol = 10^(-6); % 10^(-6) is sufficient
            x = linspace(0,1-1/N,N); % physical space domain
            x_2048 = linspace(0,1-1/2048,2048); % required for LuLu optIC interpolation to higher N
            K0 = [0:N/2-1 0 -N/2+1:-1]; % fourier space domain
            K1 = [0:N/2 -N/2+1:-1]; % fourier space domain for derivative 1
            nu = 0.001; % viscosity coefficient
            lam = [ 1, 0.5, 0.1, 0.01, 0.001];
            lambda = lam(lamtest); % smoothing parameter
            % E0 = % initial enstrophy
            % T = % length of time window
            % uField = % solution to unknown variable
            % tvector = % solution timepoints
            % testcase = [ '_' num2str(N) '_slingshot5Tfull_long' num2str(tptest) '' ]; % " _[name] " % adjust prior to run

            EnstPoints = 31;
            Enstrophy = logspace(3,6,EnstPoints);
                            
            % enspt = ensstart;
            ss_shift = 0;
            for enspt = ensstart:ensend % ss_shift = 2:3  % 
                
                %{
                if enspt > 8
                    tptest = 22;
                elseif enspt > 6
                    tptest = 21;
                elseif (enspt == 6 || enspt == 1 )
                    tptest = 19;
                elseif enspt == 5
                    tptest = 20;
                elseif (enspt == 4 || enspt == 2 )
                    tptest = 18;
                elseif enspt == 3 
                    tptest = 17;
                end
                %}

                E0 = Enstrophy(enspt);
                prefactor = 2;
                TimePoints = 50; % number of timepoints in each initial enstrophy branch
                T_ens_UB = prefactor*(1/sqrt(E0)); % use for max enstrophy T 
                T_ens_LB = T_ens_UB/TimePoints; % Adjust for decreasing T_max  
                TimeWindow = linspace(T_ens_LB,T_ens_UB,TimePoints); 
        
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
        
                for timept = timestart:timeend % 

                    ParamPoints = 2; 
                    % adjust what is being tested
                    enstest = NaN(1,ParamPoints);
                    phitest = NaN(ParamPoints,N);

                    if s == 0
                        testcase_max = [ '_' num2str(N) '_LuMax_t' num2str(stepscale) '_fine50' ];
                        testcase1 = [ '_' num2str(N) '_LuCont1_t' num2str(stepscale) '_fine50' ];
                        testcase0 = [ '_' num2str(N) '_LuCont0_t' num2str(stepscale) '_fine50' ];
                        testcase1_phi = [ '_' num2str(N) '_LuCont1_t' num2str(stepscale) '_short120' ];
                        testcase0_phi = [ '_' num2str(N) '_LuCont0_t' num2str(stepscale) '_short120' ];
                        % testcase_max = [ '_' num2str(N) '_L1Max_t' num2str(stepscale) '_tp(' num2str(tptest) ')' ];
                        % testcase1 = [ '_' num2str(N) '_L1Cont1_t' num2str(stepscale) '_tp(' num2str(tptest) ')' ]; 
                        % testcase0 = [ '_' num2str(N) '_L1Cont0_t' num2str(stepscale) '_tp(' num2str(tptest) ')' ];
                        if ss_shift ~= 0
                            testcase_max = [ '_' num2str(N) '_L1Max_t' num2str(stepscale) '_tp(' num2str(tptest) ')_shift' num2str(ss_shift) '' ];
                            testcase1 = [ '_' num2str(N) '_L1Cont1_t' num2str(stepscale) '_tp(' num2str(tptest) ')_shift' num2str(ss_shift) '' ]; 
                            testcase0 = [ '_' num2str(N) '_L1Cont0_t' num2str(stepscale) '_tp(' num2str(tptest) ')_shift' num2str(ss_shift) '' ];
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
                            CONT = 3;
                            testcase = testcase0;
                            testcase_phi = testcase0_phi;
                        end
    
                        tic
                        % initialize diagnostics
                        J = NaN(MAXITER);
                        K = NaN(MAXITER);
                        step = NaN(MAXITER);
                        momentum = NaN(MAXITER);
                        
                        % Initial condition of physical problem %
                        if ( timept == timestart || CONT ~= 1 ) 
                            CONT = 0;
                            if (CONT == 0 )%|| s == 0)
                                [ phi , ~ ] = initialguess('exact', x, E0, 0, x_2048, timept,tptest,testcase,lambda,ss_shift);
                            elseif CONT == 2
                                [ phi , ~ ] = initialguess('optIC', x, E0, 0, x_2048, timept,tptest,testcase_phi,lambda,ss_shift);
                            elseif (CONT == 3 || s > 0)
                                % shotname = ['slingshot' num2str(s) ''];
                                % shotname = 'ss_shift';
                                % shotname = 'ss_root';
                                % shotname = 'ss_start';
                                % shotname = 'slingshot';
                                [ phi , T_interval ] = initialguess(shotname, x, E0, 0, x_2048, timept,tptest,testcase,lambda,ss_shift);
                                % Adjust time interval (only when testing slingshot points)
                                if timept == timestart 
                                    % TimeWindow = linspace( T_interval/2 , 2*T_interval , timeend ); 
                                end
                            end
    
                            phi = adjust_optIC(phi,E0,K0,N);
                            phi_x = 2*pi*1i*K0.*phi;
                            E = 0.5*sum(abs(phi_x).^2)/N^2;
                        end
    
                        if ss < 7
                            T = TimeWindow(timept)*(1000+shiftT)/1000;
                        else
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
                                case 2
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
            
                        phi_save(phi, timept, E0, timestart, testcase);
                        enstest(1,p) = max(f_ens);
                        phitest(p,:) = phi;
                        runtime_save(runtime, testcase, E0, timept, lambda, timestart);
                    end %parampt

                    % choose higher enstrophy result
                    if enstest(1,1) > enstest(1,2)
                        pc = 1;
                        testcase_old = testcase1;
                        testcase_other = testcase0;
                    else
                        pc = 2;
                        testcase_old = testcase0;
                        testcase_other = testcase1;
                    end

                    % save max branch (cont or no cont)
                    phi = phitest(pc,:);
                    max_save(timept,lambda,E0,testcase_max,testcase_old,testcase_other,timeend); 

                    if timept == timeend
                        growth_save(lambda,E0,testcase_max);    
                    end

                end %timept
            end %enspt
        end %ttest
    end
return