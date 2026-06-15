function main_1DBurgers_Enstrophy(testnumber)

    % To execute, run: main_1DBurgers_Enstrophy(1);
    
    %{
    Update: Jovan Zigic, 2026/06
    Contributors: Jovan Zigic, Diego Ayala, Bartosz Protas
    Description: 
        Computing Enstrophy Growth (subject to initial enstrophy condition) via 
        Riemannian Gradient Optimization of Heuristic Choice of Initial Data for 
        1D viscous Burgers (in Fourier domain).
    Papers:
        1) Maximum Enstrophy Growth in Burgers Equation, D. Ayala, 2010 (MS Thesis)
        2) On maximum enstrophy growth in a hydrodynamic system, D. Ayala and B. Protas, 2011 (Physica D)
        3) Search for Distinct Maximizers to Finite-Time Enstrophy Growth in One-Dimensional Burgers Flows, J. Zigic, 2024 (unpublished)
    Design:
    1. Choose initial enstrophy values for computation.
        Default: E0 = 1000 <-> compute for initial enstrophy of 1000.
    2. Choose branch of time windows for E0 level.
        Default: T = [2/(50*sqrt(E0)), 2/(1*sqrt(E0))] <-> maximum enstrophy for E0 expected at T ~ 1/sqrt(E0) for E0 > 10^4.
    3. Choose continuation or no continuation of initial data for branch computation.
        p=1 uses previous timept solution, p=2 uses (Doering, Lu 2009) initial data.
    4. Search for optimal solution using gradient optimization routine.
    %}
    
    close all                                                                                             %

    root = pwd;
    addpath(genpath(fullfile(root,'data')));
    addpath(genpath(fullfile(root,'src')));
    
    % initialize diagnostic switches %
    opt_sol = 1;                            % 0: solve 1D Burgers in time, 1: optimize for maximum enstrophy                           
    Nscale = 1;                             % physical domain resolution scaling
    stepscale = 1;                          % 1/stepscale = time-step scaling for oscillations in enstrophy values
    VERBOSE = 1;                            % set level of data output
    
    % declare and initialize parameters %
    ascent = 'RCGPR';                       % Riemmanian Conjugate Gradient with Polak-Ribiere
    CONTSET = 0;                            % ( 0: Lu, 2: other data file)
    EnstPoints = 1;                        % number of initial enstrophy points
    Enstrophy = logspace(2,3,EnstPoints);   %
    N = 2048*Nscale;                        % number of grid points
    MAXITER = 1000;                         % 1000; % this should never be reached
    J_step_tol = 10^(-6);                   % 10^(-6) is sufficient
    x = linspace(0,1-1/N,N);                % physical space domain
    x_2048 = linspace(0,1-1/2048,2048);     % required for LuLu optIC interpolation to higher N
    K0 = [0:N/2-1 0 -N/2+1:-1];             % fourier space domain
    K1 = [0:N/2 -N/2+1:-1];                 % fourier space domain for derivative 1
    nu = 0.001;                             % viscosity coefficient
    lambda = 1;                             % smoothing parameter  
    phi = 0;                                % optimal initial data
    
    %{
    E0      : initial enstrophy
    T       : length of time window
    uField  : solution to unknown variable
    tvector : solution timepoints
    %}
    
    for enspt = 1:length(Enstrophy) % "enspt" <-> enstrophy point
    
        % choose time windows to test at each E0 %
        E0 = Enstrophy(enspt);                                  %
        prefactor = 2;                                          %
        TimePoints = 3;                                        % number of timepoints in each initial enstrophy branch
        T_ens_UB = prefactor*(1/sqrt(E0));                      % use for max enstrophy T 
        T_ens_LB = T_ens_UB/TimePoints;                         % Adjust for decreasing T_max  
        TimeWindow = linspace(T_ens_LB,T_ens_UB,TimePoints);    %
    
        % choose data output at each E0 %
        switch VERBOSE                                                          %
            case 1                                                              %
                if not(isfolder([pwd  '/data/diagnostics' ]))                       % create local directories for data storage
                    mkdir([pwd  '/data/diagnostics' ]);                             %
                    mkdir([pwd  '/data/enstrophy_solution' ]);                      %
                    mkdir([pwd  '/data/time_evolution' ]);                          %
                    mkdir([pwd  '/data/spectrum' ]);                                %
                    mkdir([pwd  '/data/spectrum/temp' ]);                           %
                end
            case 3                                                          %
                if not(isfolder([pwd  '/data/kappa' ]))                         % create local directories for data storage
                    mkdir([pwd  '/data/kappa' ]);                               %
                    mkdir([pwd  '/data/kappa/kappa_E0_' num2str(E0) '' ]);      %
                end  
        end

        if not(isfolder([pwd  '/data/spectrum/spectrum_E0_' num2str(E0) '' ]))                                           % create local directories for data storage
                    mkdir([pwd  '/data/spectrum/spectrum_E0_' num2str(E0) '' ]);    %
        end
    
        for timept = 1:length(TimeWindow)          %
    
            T = TimeWindow(timept);
            if opt_sol == 1
                ParamPoints = 2;                % test both continuation and data file
            else
                ParamPoints = 1;                % test data file
            end
            enstest = NaN(1,ParamPoints);       %
            phitest = NaN(ParamPoints,N);       %
    
            %%% name testcase for output data %%%
            [testcase_max,testcase1,testcase0] = testcase_name(N,stepscale,ascent,testnumber);
    
            for p = 1:ParamPoints         
        
                if p == 1 && opt_sol == 1   %
                    CONT = 1;               %
                    testcase = testcase1;   %
                else
                    CONT = CONTSET;         %
                    testcase = testcase0;   %
                end
    
                tic                         % start timer
                
                % solve optimization problem if opt_sol = 1 (otherwise,
                % solve forward problem)
                [phi,f_ens] = optimize_Burgers(opt_sol,CONT,testcase,enspt,timept,phi,N,T,...
                                                x,Enstrophy,TimeWindow,stepscale,ascent,VERBOSE,...
                                                E0,x_2048,lambda,J_step_tol,K0,K1,nu,MAXITER);

                phi_save(phi, timept, E0, testcase);         % save optimal initial data 
                fprintf('Saved optimized initial data file at %01dh%02dm%02ds\t \n',floor(toc/3600),floor(mod(toc/60,60)),floor(mod(toc,60))')%
                enstest(1,p) = max(f_ens);                              % save max enstrophy value
                phitest(p,:) = phi;                                     % save optimal initial data for continuation
                %runtime_save(runtime, testcase, E0, timept, lambda);    % save runtime data for each routine
    
            end  
            
            % choose higher enstrophy result %
            if opt_sol == 1
                if enstest(1,1) > enstest(1,2) || isnan(enstest(1,2))   %
                    pc = 1;                                             %
                    testcase_old = testcase1;                           %
                    testcase_other = testcase0;                         %
                else                                %
                    pc = 2;                         %
                    testcase_old = testcase0;       %
                    testcase_other = testcase1;     %
                end
        
                % save data with higher enstrophy (continuation or no continuation) %
                phi = phitest(pc,:);                                                            %
                max_save(timept,lambda,E0,testcase_max,testcase_old,testcase_other,TimeWindow);    %
    
                if timept == length(TimeWindow)                        %
                    growth_save(lambda,E0,testcase_max);    % save max enstrophy across branch of time windows   
                end
            end
    
        end %timept
    end %enspt
return