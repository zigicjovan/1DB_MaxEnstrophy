function main_1DBEnstrophy_2024_v4_1(ascent,caseID)

%{
[v4_1: update diagnostics metrics, adjust stopping criteria, remove non-Riemannian methods]
Update: Jovan Zigic, 2024/06
Contributors: Jovan Zigic, Diego Ayala, Bartosz Protas
Description: 
    Computing Enstrophy Growth (subject to initial enstrophy condition) via 
    Riemannian Gradient Optimization of Heuristic Choice of Initial Data for 
    1D viscous Burgers (in Fourier domain).
Papers:
    1) Maximum Enstrophy Growth in Burgers Equation, D. Ayala, 2010 (MS Thesis)
    2) On maximum enstrophy growth in a hydrodynamic system, D. Ayala and B. Protas, 2011 (Physica D)
    3) Search for Distinct Maximizers to Finite-Time Enstrophy Growth in
    One-Dimensional Burgers Flows, J. Zigic, 2024 (unpublished)
    4) Momentum-Based Gradient Descent Methods for PDE-Constrained
    Optimization on Riemannian Manifolds, J. Zigic and B. Protas, 2024 (in
    progress)
Design (by level) - v3 update:
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
5. Choose continuation or no continuation of initial data for branch computation and compute in parallel.
    Default: p = 1:2 <-> p=1 uses previous timept solution, p=2 uses (Doering, Lu 2009) initial data.
6. Search for optimal solution using gradient optimization routine
%}

    s = 0; % slingshot number (0 if not using slingshot initial data)
    shiftsign = 0; % choose whether to shift time window (0: no shift, 1: positive shift, 2: negative shift)
    lamtest = 3; % choice of smoothing parameter

    for liptest = 1:7 % choice of smoothing parameter

        % sensitivity testing of time window perturbation % 
        if shiftsign < 2
            pm = 'p'; % positive shift
        else
            pm = 'm';  % negative shift
        end

        % shiftT = 0; % time window shift parameter
        tptest = 1; % slingshot search parameter

        for shiftT = 0:0 % "shiftT" <-> sensitivity testing || "tptest" <-> slingshot searching
            clearvars -except shiftT tptest s shiftsign pm lamtest ascent caseID liptest condtest; 
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
            timeend = 10;% 

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
            lipparameter = [ 150, 170, 190, 200, 210, 230, 250 ];
            Lip = lipparameter(liptest);
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
                TimeWindow = linspace(T_ens_LB,T_ens_UB,TimePoints); 
        
                % choose data output at each E0 %
                switch VERBOSE
                    case 1 
                        mkdir([pwd  '/data/diagnostics' ]);
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
                    psitest = NaN(ParamPoints,N);
                    Etest = NaN(1,ParamPoints);

                    %%% name testcase for output data %%%
                    [testcase_max,testcase1,testcase0] = testcase_name(caseID,N,stepscale,ascent,tptest,pm,s,shiftT,ss_shift,0);

                    for p = 1:ParamPoints % do in parallel
                
                        if p == 1
                            CONT = 1;
                            testcase = testcase1;
                        elseif p == 2
                            CONT = CONTSET;
                            testcase = testcase0;
                        end

                        if timept > timestart
                            phi = load('p_vars.mat','phi').phi;
                            psi = load('p_vars.mat','psi_save').psi_save;
                            E = load('p_vars.mat','E_save').E_save;
                        end
    
                        tic % start timer
                        % initialize diagnostics
                        J = NaN(MAXITER,1);
                        Jeval = NaN(5*MAXITER,1);
                        K = NaN(MAXITER,1);
                        step = NaN(MAXITER,1);
                        momentum = NaN(MAXITER,1);
                        
                        % Initial condition of physical problem %
                        if ( timept == timestart || CONT ~= 1 ) 
                            if ( CONT == 0 || CONT == 1 )
                                [ phi , ~ ] = initialguess('exact', x, E0, 0, x_2048, timept,tptest,testcase,lambda,ss_shift);
                                % [ phi , ~ ] = initialguess('sine', x, E0, 0, x_2048, timept,tptest,testcase,lambda,ss_shift);
                                % [ phi , ~ ] = initialguess('sin2', x, E0, 0, x_2048, timept,tptest,testcase,lambda,ss_shift);
                                % [ phi , ~ ] = initialguess('2 modes', x, E0, 0, x_2048, timept,tptest,testcase,lambda,ss_shift);
                            elseif CONT == 2
                                [ phi , ~ ] = initialguess('optIC', x, E0, 0, x_2048, timept,tptest,testcase,lambda,ss_shift);
                            elseif CONT == 3 
                                %%% choose slingshot setting %%%
                                % shotname = ['slingshot' num2str(s) ''];
                                shotname = 'ss_shift';
                                % shotname = 'ss_root';
                                % shotname = 'ss_start';
                                % shotname = 'slingshot';
                                [ phi , ~ ] = initialguess(shotname, x, E0, 0, x_2048, timept,tptest,testcase,lambda,ss_shift);
                                % adjust time interval only when testing slingshot points
                                if timept == timestart
                                    % TimeWindow = linspace( T_interval/2 , 2*T_interval , timeend ); 
                                end
                            end
    
                            phi = retraction(phi,E0,K0,N); % initial data 
                            phi_x = 2*pi*1i*K0.*phi; % derivative of initial data
                            E = 0.5*sum(abs(phi_x).^2)/N^2; % enstrophy E0 for retraction
                            psi = phi; % projection term
                        end
    
                        if shiftsign < 2
                            T = TimeWindow(timept)*(1000+shiftT)/1000;
                        else
                            T = TimeWindow(timept)*(1000-shiftT)/1000;
                        end

                        ITER = 1; % start counting iterations
                        iterJ = 0; % start counting iterations of function evaluation
                        
                        fprintf('\n\n  Enstrophy point = %d,  Time point = %d \n',enspt,timept);
                        fprintf('  Iterations: ')
    
                        % Do at least one iteration
                        [tvector,uField] = BurgersDS_Fourier(phi,K1,K0,T,nu,N,stepscale); % solve PDE
                        Ntime = length(tvector); % number of time-steps in PDE solution 
                        uu = uField(Ntime,:); % physical solution at final time-step
                        J(ITER) = eval_J(uu,phi,K0,N); % evaluate objective functional (E(T) - E0)
                        iterJ = iterJ + 1;
                        Jeval(iterJ) = J(ITER);
                        K(ITER) = norm(uu); % kinetic energy at E(T)
                        step(ITER) = 0;
                        momentum(ITER) = 0;
                        
                        u_adj = BurgersAS_Fourier(uu,K1,K0,T,nu,N,uField,tvector,stepscale); % solve adjoint PDE
                        gradJ = eval_grad_J(u_adj,phi,K1,lambda); % gradient of objective functional
                        gradJ_x = (2*pi*1i*K0).*gradJ; % derivative of gradient of objective functional
                        phi_x = 2*pi*1i*K0.*phi; % derivative of physical data
                        
                        %%% Gradient optimization initial step %%%
                        switch caseID
                            case 'p40'
                                cond = 0.40; % conditioning of HB/NAG
                            case 'p42'
                                cond = 0.42; % conditioning of HB/NAG
                            case 'p44'
                                cond = 0.44; % conditioning of HB/NAG
                            case 'p46'
                                cond = 0.46; % conditioning of HB/NAG
                            case 'p48'
                                cond = 0.48; % conditioning of HB/NAG
                            case 'p50'
                                cond = 0.50; % conditioning of HB/NAG
                        end
                        % cond = 0.4; % conditioning of HB/NAG
                        % Lip = 200; % Lipschitz constant approximation ~ upper bound of quadratic growth
                        mu_bound = cond*Lip; % approximation ~ lower bound of quadratic growth
                        cond_inv = mu_bound/Lip; % inverse of condition number
                        switch ascent
                            case {'RCGPR','RCGRMIL'}
                                %%% Riemannian Conjugate Gradient algorithm %%%
                                tau = -2*real( sum(phi.*conj(gradJ)) + lambda^2*sum(phi_x.*conj(gradJ_x)) ) ... 
                                    /( sum(abs(gradJ).^2) + lambda^2*sum(abs(gradJ_x).^2) ) ; % initial guess for optimal step length
                                proj_gradJ = projection(gradJ,psi,K0,N,lambda); % projected gradient
                                direction = proj_gradJ; % direction towards optimal solution
                                [tau,Jeval,iterJ] = optimize_tau(phi,direction,abs(0.2*tau),E,K1,K0,T,nu,N,stepscale,Jeval,iterJ); % compute initial optimal step size 
                                phi = phi + tau*direction; % iterative solution update
                                psi = phi; % projection term
                            case 'RHB'
                                %%% Riemannian Heavy-Ball algorithm %%%
                                alpha = ( 2 / ( sqrt(Lip) + sqrt(mu_bound) ) )^2;
                                beta = ( ( 1 - sqrt(cond_inv) ) / ( 1 + sqrt(cond_inv) ) )^2;
                                proj_gradJ = projection(gradJ,psi,K0,N,lambda); % projected gradient
                                direction = alpha*proj_gradJ; % direction towards optimal solution
                                phi_old = phi;
                                phi = phi + direction; % iterative solution update
                                psi = phi; % projection term
                            case 'RNAG'
                                %%% Riemannian Nesterov Accelerated Gradient algorithm %%%
                                alpha = 4 / ( 3*Lip + mu_bound) ; % Trung Vu (2018) coefficients
                                beta = ( 2 - sqrt(3*cond_inv + 1) ) / ( 2 + sqrt(3*cond_inv + 1) ) ; % Trung Vu (2018) coefficients
                                proj_gradJ_int = projection(gradJ,psi,K0,N,lambda); % projected gradient
                                direction = alpha*proj_gradJ_int; % direction towards optimal solution
                                phi_old = phi;
                                phi = phi + direction; % iterative solution update
                                psi = phi; % projection term
                        end
                        phi = retraction(phi,E,K0,N); % retraction to constraint manifold
                        %%%
                        
                        [tvector,uField] = BurgersDS_Fourier(phi,K1,K0,T,nu,N,stepscale);
                        Ntime = length(tvector);
                        uu = uField(Ntime,:);
                        J_new = eval_J(uu, phi,K0,N);
                        iterJ = iterJ + 1;
                        Jeval(iterJ) = J_new;
                        
                        J_step = (J_new - J(ITER))/J(ITER); % relative change in objective functional value
                        ITER = ITER + 1;
                        J(ITER) = J_new; 
                        K(ITER) = norm(uu);
                        switch ascent
                                case {'RCGPR','RCGRMIL'}
                                    step(ITER) = tau;
                                case {'RHB','RNAG'}
                                    step(ITER) = alpha;
                        end
                        momentum(ITER) = 0;
        
                        switch VERBOSE
                            case 3  % test definition of Gateaux differential
                                u_adj = BurgersAS_Fourier(uu,K1,K0,T,nu,N,uField,tvector,stepscale);
                                gradJ_new = eval_grad_J(u_adj,phi,K1,lambda);
                                gradJ_x_new = (2*pi*1i*K0).*gradJ_new;
                                gradJ = gradJ_new;
                                gradJ_x = gradJ_x_new;
                                [epsilon,kappa] = kappaTestFourier(phi,gradJ,gradJ_x,J(ITER),x,K1,K0,T,nu,N,lambda,stepscale);
                                kappa_save(kappa, epsilon, timept,E0,lambda,ITER,testcase);
                        end
                    
                        while ( J_step > J_step_tol ) && ( ITER <= MAXITER )
                            
                            %%% gradient momentum adjustment %%%
                            switch ascent
                                case {'RCGPR'}
                                    beta_denom = sum(abs(gradJ).^2)/N^2 + lambda^2*sum(abs(gradJ_x).^2)/N^2; % Polak-Ribiere denominator
                                case {'RCGRMIL'}
                                    dir_x = (2*pi*1i*K0).*direction; % derivative of conjugate direction
                                    beta_denom = sum(abs(direction).^2)/N^2 + lambda^2*sum(abs(dir_x).^2)/N^2; % RMIL denominator
                            end
                            u_adj = BurgersAS_Fourier(uu,K1,K0,T,nu,N,uField,tvector,stepscale);
                            gradJ_new = eval_grad_J(u_adj,phi,K1,lambda);
                            gradJ_x_new = (2*pi*1i*K0).*gradJ_new;
                            beta_num = sum(gradJ_new.*conj(gradJ_new - gradJ))/N^2 ...
                                + lambda^2*sum(gradJ_x_new.*conj(gradJ_x_new - gradJ_x))/N^2; % PR/RMIL numerator
                            gradJ = gradJ_new;
                            gradJ_x = gradJ_x_new;
                            switch ascent
                                case {'RCGPR','RCGRMIL'}
                                    % reset beta after 20 iterations %
                                    if mod(ITER,20)==0
                                        beta=0;
                                    else
                                        beta = max(real(beta_num/beta_denom), 0);
                                    end
                            end
                            %%%
                            
                            %%% Gradient optimization iterative step %%%
                            switch ascent
                                case {'RCGPR','RCGRMIL'}
                                    proj_gradJ = projection(gradJ,psi,K0,N,lambda);
                                    psi_x = (2*pi*1i*K0).*psi; % derivative of projector term
                                    trans_dir = projection(direction,psi,K0,N,lambda) ...
                                        /( sum(abs(psi).^2) + lambda^2*sum(abs(psi_x).^2) ); % vector transport of previous direction
                                    direction = proj_gradJ + beta*trans_dir; % direction towards optimal solution
                                    [tau,Jeval,iterJ] = optimize_tau(phi,direction,tau,E,K1,K0,T,nu,N,stepscale,Jeval,iterJ); % compute optimal step size 
                                    phi = phi + tau*direction; % iterative solution update
                                    psi = phi; % projection term
                                case 'RHB'
                                    proj_gradJ = projection(gradJ,psi,K0,N,lambda);
                                    trans_dir = projection(direction,psi,K0,N,lambda); % projected previous direction
                                    direction = alpha*proj_gradJ + beta*trans_dir; % direction towards optimal solution
                                    phi_old = phi;
                                    phi = phi + direction; % iterative solution update
                                    psi = phi; % projection term
                                case 'RNAG'
                                    d_int = phi - phi_old;
                                    trans_dir = projection(d_int,psi,K0,N,lambda); % projected previous direction
                                    phi_int = phi + beta*trans_dir;
                                    phi_int = retraction(phi_int,E,K0,N);
                                    psi = phi_int; % new projection term
                                    [tvector_int,uField_int] = BurgersDS_Fourier(phi_int,K1,K0,T,nu,N,stepscale);
                                    Ntime_int = length(tvector_int);
                                    uu_int = uField_int(Ntime_int,:);
                                    u_adj_int = BurgersAS_Fourier(uu_int,K1,K0,T,nu,N,uField_int,tvector_int,stepscale);
                                    gradJ_int = eval_grad_J(u_adj_int,phi_int,K1,lambda);
                                    proj_gradJ_int = projection(gradJ_int,psi,K0,N,lambda); % projected gradient
                                    direction = alpha*proj_gradJ_int; % direction towards optimal solution
                                    phi_old = phi;
                                    phi = phi + direction; % iterative solution update
                                    psi = phi; % projection term
                            end
                            phi = retraction(phi,E,K0,N); % retraction to constraint manifold
                            %%%
                            
                            [tvector,uField] = BurgersDS_Fourier(phi,K1,K0,T,nu,N,stepscale);
                            Ntime = length(tvector);
                            uu = uField(Ntime,:);
                            J_new = eval_J(uu, phi,K0,N);
                            iterJ = iterJ + 1;
                            Jeval(iterJ) = J_new;
                            
                            J_step = (J_new - J(ITER))/J(ITER);
                            ITER = ITER+1;
                            J(ITER) = J_new;
                            K(ITER) = norm(uu);
                            switch ascent
                                case {'RCGPR','RCGRMIL'}
                                    step(ITER) = tau;
                                case {'RHB','RNAG'}
                                    step(ITER) = alpha;
                            end
                            momentum(ITER) = beta;
        
                            disp(['P = ' num2str(p) ' Completed iteration ' num2str(ITER-1)]);
                            disp(['P = ' num2str(p) ' Current runtime: ' num2str(toc)]);
                        end
        
                        if J_step < 0
                            phi = phi_old; % cancel out negative step
                        end

                        runtime = toc; % stop timer
                        
                        disp('');
                        disp(['P = ' num2str(p) ' Time window: ' num2str(T)]);
                        disp(['P = ' num2str(p) ' Initial Enstrophy: ' num2str(E0)]);
                        disp(['P = ' num2str(p) ' Number of iterations: ' num2str(ITER-1)]);
                        disp(['P = ' num2str(p) ' Total Runtime: ' num2str(runtime)]);
                        disp(['P = ' num2str(p) ' Final Relative Step: ' num2str(J_step)]);
                        
                        if ITER > MAXITER
                            disp('Optimization reached maximum number of iterations.');
                            disp(['P = ' num2str(p) ' Final Relative Step: ' num2str(J_step)]);
                        end
            
                        templam = lambda; % keep this until Lip and mu_bound are chosen
                        lambda = Lip; % keep this until Lip and mu_bound are chosen
                        switch VERBOSE
                            case 1 % save data
                                diagnostics_save(J,K,step,momentum,timept,lambda,E0,ITER,testcase);
                                evaluation_save(Jeval,timept,lambda,E0,iterJ,testcase);
                                f_ens = time_evolution_and_optIC_save(Ntime,phi,timept,E0,lambda,T,N,uField,testcase,K0,stepscale);
                                branch_save(timept,f_ens,E0,lambda,T,testcase); 
                        end
            
                        phi_save(phi, timept, E0, timestart, testcase); % save optimal initial data 
                        enstest(1,p) = max(f_ens); % save max enstrophy value
                        phitest(p,:) = phi; % save optimal initial data for continuation
                        psitest(p,:) = psi;
                        Etest(1,p) = E;
                        runtime_save(runtime, testcase, E0, timept, lambda); % save runtime data for each routine

                        lambda = templam; % keep this until Lip and mu_bound are chosen

                    end %parampt
                    % end parallel

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

                    templam = lambda;
                    lambda = Lip;
                    % save data with higher enstrophy (continuation or no continuation) %
                    phi = phitest(pc,:);
                    max_save(timept,lambda,E0,testcase_max,testcase_old,testcase_other,timeend); 

                    if timept == timeend
                        growth_save(lambda,E0,testcase_max); % save max enstrophy across branch of time windows   
                    end

                    lambda = templam;

                    psi_save = psitest(pc,:);
                    E_save = Etest(1,pc);

                    save('p_vars.mat','phi','psi_save','E_save');

                end %timept
            end %enspt
        end 
    end
return