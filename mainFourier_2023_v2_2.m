function mainFourier_2023_v2_2

% Version control %
%{
v1_0 : added VERBOSE and master diagnostics
v1_1 : time evolution, branch and spectrum
v1_2 : spectrum evolution, E0 starts at 10e3
v1_3 : T ~ C*(1/sqrt(E0)), higher resolution, plot E0 vs {Emax,T} growth
v2_0 : cleanup of code, change TimeWindow bounds, kappa test
v2_1 : fixed growth issues, global cleanup, increase tolerance
v2_2 : add runtime, manual time_start, IC correction (in adjust_optIC)
%}

    clearvars; close all;

    % initialize diagnostic switches %
    VERBOSE = 1; 
    testcase = '_6144_optIC2048ContLu7_t5b'; % " _[name] " % adjust prior to run

    % declare and initialize parameters %
    CONT = 2; % no continuation recommended (0 for 2048 opt IC, 1 for cont, 2 for Lu)
    stepscale = 5; % 1/stepscale = time-step scaling for oscillations in enstrophy values
    MAXITER = 6; % 6 iterations is sufficient
    J_step_tol = 10^(-6); % 10^(-6) is sufficient
    N = 2048*3; % number of grid points
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
    
    ParamPoints = 2;
      
    for p = 2:ParamPoints
    
        % number of timepoints in each initial enstrophy branch
        T_res = 1; 
        TimePoints = 31*ceil(T_res);
        % amplitude = linspace(0,10/pi,5);

        if p == 1
            EnstPoints = 31;
            Enstrophy = logspace(0,3,EnstPoints);
            EnstPoints = 30;
            ens_start = 21;
        elseif p == 2
            EnstPoints = 31;
            Enstrophy = logspace(3,6,EnstPoints);
            EnstPoints = 16;
            ens_start = 16;
        end

        time_start = 13;
        
        for enspt = ens_start:EnstPoints
            
            % CONT = 0;
            E0 = Enstrophy(enspt);
            prefactor = 2;
            T_ens_UB = prefactor*(1/sqrt(E0)); % Adjust for decreasing T_max    
            T_ens_LB = T_ens_UB/TimePoints; % Adjust for decreasing T_max  
            TimeWindow = linspace(T_ens_LB,T_ens_UB,TimePoints); 
    
            switch VERBOSE
                case 1 
                    mkdir([pwd  '/data/diagnostics/diagnostics_E0_' num2str(E0) '' ]);
                    mkdir([pwd  '/data/enstrophy_solution' ]);
                    mkdir([pwd  '/data/time_evolution' ]);
                    mkdir([pwd  '/data/spectrum' ]);
                    mkdir([pwd  '/data/spectrum/spectrum_E0_' num2str(E0) '' ]);
                    mkdir([pwd  '/data/kappa/kappa_E0_' num2str(E0) '' ]);
                    mkdir([pwd  '/data/runtime' ]);
                case 2
                    mkdir([pwd  '/data/enstrophy_solution' ]);
                    mkdir([pwd  '/data/time_evolution' ]);
                    mkdir([pwd  '/data/spectrum' ]);
                    mkdir([pwd  '/data/spectrum/spectrum_E0_' num2str(E0) '' ]);
                    mkdir([pwd  '/data/runtime' ]);
                case 3
                    mkdir([pwd  '/data/kappa' ]);
            end
    
            for timept = time_start:TimePoints % 

                %{
                if enspt == 16 %32
                    if timept == 13
                        CONT = 1;
                    elseif timept > 14
                        CONT = 1;
                    else 
                        CONT = 0;
                    end
                elseif enspt == 17 %38
                    if timept > 17
                        CONT = 1;
                    elseif timept == 17
                        CONT = 0;
                    elseif timept > 12 
                        CONT = 1;
                    end
                elseif enspt == 18 %50
                    if timept > 12
                        CONT = 1;
                    end
                elseif enspt == 19 %63
                    if timept > 14
                        CONT = 1;
                    elseif timept == 14
                        CONT = 0;
                    elseif timept > 10 
                        CONT = 1;
                    end
                elseif enspt == 20 %79
                    if timept > 15
                        CONT = 1;
                    elseif timept == 15
                        CONT = 0;
                    elseif timept == 14 
                        CONT = 1;
                    elseif timept == 13
                        CONT = 0;
                    elseif timept > 10 
                        CONT = 1;
                    end
                elseif enspt == 21 %100
                    if timept > 16
                        CONT = 1;
                    elseif timept > 14
                        CONT = 0;
                    elseif timept == 14 
                        CONT = 1;
                    elseif timept == 13
                        CONT = 0;
                    elseif timept > 10 
                        CONT = 1;
                    elseif timept == 10 
                        CONT = 0;
                    elseif timept == 9 
                        CONT = 1;
                    end
                else
                    CONT = 0;
                end
                %}

                tic
                % initialize diagnostics
                J = NaN(MAXITER);
                K = NaN(MAXITER);
                step = NaN(MAXITER);
                momentum = NaN(MAXITER);

                T = TimeWindow(timept);
                ITER = 1;
                
                fprintf('\n\n  Enstrophy point = %d,  Time point = %d \n',enspt,timept);
                fprintf('  Iterations: ')
                
                % Initial condition of physical problem %
                if ( timept == time_start || CONT ~= 1 )    
                    if ( N == 2048 || CONT == 2 ) 
                        phi = initialguess('exact', x, E0, 0, x_2048, timept,time_start);
                    else
                        phi = initialguess('optIC', x, E0, 0, x_2048, timept,time_start);
                    end
                end
                
                phi = adjust_optIC(phi,E0,K0,N);
                phi_x = 2*pi*1i*K0.*phi;
                E = 0.5*sum(abs(phi_x).^2)/N^2;
                
                % Do at least one iteration
                [tvector,uField] = BurgersDS_Fourier(phi,K1,K0,T,nu,N,stepscale);
                Ntime = length(tvector);
                uu = uField(Ntime,:);
                J(ITER) = eval_J(uu,phi,K0,N); % iter case 1
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
                    disp(['Runtime so far: ' num2str(toc)]);
                end

                runtime = toc;
                
                disp('');
                disp(['Lambda: ' num2str(lambda)]);
                disp(['Time window: ' num2str(T)]);
                disp(['Initial Enstrophy: ' num2str(E0)]);
                disp(['Number of iterations: ' num2str(ITER-1)]);
                disp(['Runtime: ' num2str(runtime)]);
                
                if ITER > MAXITER
                    disp('Simulation reached maximum number of iterations.');
                    disp(['Final Relative Objective Step: ' num2str(J_step)]);
                end
    
                switch VERBOSE
                    case 1 
                        diagnostics_save(J,K,step,momentum,timept,lambda,E0,ITER,testcase);
                        f_ens = time_evolution_and_optIC_save(Ntime,phi,timept,E0,lambda,T,N,uField,testcase,K0,time_start,stepscale);
                        branch_save(timept,f_ens,E0,lambda,T,testcase); 
                    case 2 
                        f_ens = time_evolution_and_optIC_save(Ntime,phi,timept,E0,lambda,T,N,uField,testcase,K0,time_start,stepscale);
                        branch_save(timept,f_ens,E0,lambda,T,testcase); 
                    case 3 
                        %
                end
    
                phi_save(phi, timept, E0, time_start, testcase);
                runtime_save(runtime, testcase, E0, timept, lambda, time_start);
            end
            growth_save(lambda,E0,testcase);  
        end
    end
return