function [phi,f_ens] = optimize_Burgers(opt_sol,CONT,testcase,enspt,timept,phi,N,T,...
                                        x,Enstrophy,TimeWindow,stepscale,ascent,VERBOSE,...
                                        E0,x_2048,lambda,J_step_tol,K0,K1,nu,MAXITER)

    J = NaN(MAXITER,1);         %
    Jeval = NaN(5*MAXITER,1);   %
    K = NaN(MAXITER,1);         %
    step = NaN(MAXITER,1);      %
    momentum = NaN(MAXITER,1);  %
    
    % Initial condition of physical problem %
    if ( timept == 1 || CONT ~= 1 )     %
        if ( CONT == 0 || CONT == 1 )           %
            phi = initialguess('exact', x, E0, 0, x_2048, timept,testcase,lambda);
        elseif CONT == 2                        %
            phi = initialguess('optIC', x, E0, 0, x_2048, timept,testcase,lambda);
        end

        phi = retraction(phi,E0,K0,N);          % initial data 
    end    
    psi = phi;                                  % projection term
    phi_x = 2*pi*1i*K0.*phi;                % derivative of initial data
    E = 0.5*sum(abs(phi_x).^2)/N^2;         % enstrophy E0 for retraction 

    ITER = 1;           % start counting iterations
    iterJ = 0;          % start counting iterations of function evaluation
    
    fprintf('\n Solving 1D Burgers, N = %d, E0 = %.0f (point %d of %d), T = %.4f (point %d of %d) \n', ...
        N , E0 , enspt , length(Enstrophy) , T , timept , length(TimeWindow) );        %

    % Solve 1D Burgers forward in time
    [tvector,uField] = BurgersDS_Fourier(phi,K1,K0,T,nu,N,stepscale);       % solve PDE
    Ntime = length(tvector);                                                % number of time-steps in PDE solution 
    uu = uField(Ntime,:);                                                   % physical solution at final time-step
    J(ITER) = eval_J(uu,phi,K0,N);                                          % evaluate objective functional (E(T) - E0)
    iterJ = iterJ + 1;                                                      %
    Jeval(iterJ) = J(ITER);                                                 %
    K(ITER) = norm(uu);                                                     % kinetic energy at E(T)
    step(ITER) = 0;                                                         %
    momentum(ITER) = 0;                                                     %

    if opt_sol == 1

        % Do at least one iteration
        fprintf('Iter \t J0 \t Adjoint solved \t Line-search \t Forward solved \t J1 \t\t J change \n ')
        fprintf('----------------------------------------------------------------------------------------------------------\n')
    
        fprintf('%02d\t%.4f\t', ITER, J(ITER))
        
        u_adj = BurgersAS_Fourier(uu,K1,K0,T,nu,N,uField,tvector,stepscale);    % solve adjoint PDE
        gradJ = eval_grad_J(u_adj,phi,K1,lambda);                               % gradient of objective functional
        gradJ_x = (2*pi*1i*K0).*gradJ;                                          % derivative of gradient of objective functional
        phi_x = 2*pi*1i*K0.*phi;                                                % derivative of physical data
    
        fprintf('%01dh%02dm%02ds\t',floor(toc/3600),floor(mod(toc/60,60)),floor(mod(toc,60)))
    
        %%% Save optimization steps
        save_gradJ = NaN(1,N);                      %
        save_gradJ(ITER,1:length(gradJ)) = gradJ;   %
        save_direction = zeros(1,N);                %
        
        %%% Gradient optimization initial step %%%
        switch ascent                                                                                                   %
            case {'RCGPR'}     
                %%% Riemannian Conjugate Gradient algorithm %%%
                tau = -2*real( sum(phi.*conj(gradJ)) + lambda^2*sum(phi_x.*conj(gradJ_x)) ) ... 
                    /( sum(abs(gradJ).^2) + lambda^2*sum(abs(gradJ_x).^2) ) ;                                           % initial guess for optimal step length
                proj_gradJ = projection(gradJ,psi,K0,N,lambda);                                                         % projected gradient
                direction = proj_gradJ;                                                                                 % direction towards optimal solution
                [tau,Jeval,iterJ] = optimize_tau(phi,direction,abs(0.2*tau),E,K1,K0,T,nu,N,stepscale,Jeval,iterJ);      % compute initial optimal step size 
                phi = phi + tau*direction;                                                                              % iterative solution update
                phi = retraction(phi,E,K0,N);                                                                           % retraction to constraint manifold
                psi = phi;                                                                                              % projection term
                fprintf('%01d, %02dh%02dm%02ds\t\t',iterJ,floor(toc/3600),floor(mod(toc/60,60)),floor(mod(toc,60)))
        end                             %
        
        [tvector,uField] = BurgersDS_Fourier(phi,K1,K0,T,nu,N,stepscale);   %
        Ntime = length(tvector);                                            %
        uu = uField(Ntime,:);                                               %
        J_new = eval_J(uu, phi,K0,N);                                       %
        iterJ = iterJ + 1;                                                  %
        Jeval(iterJ) = J_new;                                               %
    
        fprintf('%01dh%02dm%02ds\t',floor(toc/3600),floor(mod(toc/60,60)),floor(mod(toc,60)))
        
        J_step = (J_new - J(ITER))/J(ITER);         % relative change in objective functional value
        ITER = ITER + 1;                            %
        J(ITER) = J_new;                            %
        K(ITER) = norm(uu);                         %
        switch ascent                                                                                                       %
            case {'RCGPR'}     %
                step(ITER) = tau;                                                                                       %                                                                                   %
        end                     %
        momentum(ITER) = 0;     %
    
        fprintf('%.4f\t%d\n%02d\t%.4f\t', J(ITER), abs(J_step), ITER, J(ITER))
    
        switch VERBOSE                                                                                          %
            case 3                                                                                              % test definition of Gateaux differential
                u_adj = BurgersAS_Fourier(uu,K1,K0,T,nu,N,uField,tvector,stepscale);                            %
                gradJ_new = eval_grad_J(u_adj,phi,K1,lambda);                                                   %
                gradJ_x_new = (2*pi*1i*K0).*gradJ_new;                                                          %
                gradJ = gradJ_new;                                                                              %
                gradJ_x = gradJ_x_new;                                                                          %
                [epsilon,kappa] = kappaTestFourier(phi,gradJ,gradJ_x,J(ITER),x,K1,K0,T,nu,N,lambda,stepscale);  %
                kappa_save(kappa, epsilon, timept,E0,lambda,ITER,testcase);                                     %
        end
    
        while ( abs(J_step) > J_step_tol ) && ( ITER <= MAXITER )           %
            
            %%% gradient momentum adjustment %%%
            switch ascent                                                                                       %
                case {'RCGPR'}   %
                    beta_denom = sum(abs(gradJ).^2)/N^2 + lambda^2*sum(abs(gradJ_x).^2)/N^2;                    % Polak-Ribiere denominator
            end
            u_adj = BurgersAS_Fourier(uu,K1,K0,T,nu,N,uField,tvector,stepscale);    %
            gradJ_new = eval_grad_J(u_adj,phi,K1,lambda);                           %
            gradJ_x_new = (2*pi*1i*K0).*gradJ_new;                                  %
            beta_num = sum(gradJ_new.*conj(gradJ_new - gradJ))/N^2 ...              %
                + lambda^2*sum(gradJ_x_new.*conj(gradJ_x_new - gradJ_x))/N^2;       % PR numerator
            gradJ = gradJ_new;                                                      %
            gradJ_x = gradJ_x_new;                                                  %
    
            fprintf('%01dh%02dm%02ds\t',floor(toc/3600),floor(mod(toc/60,60)),floor(mod(toc,60)))
    
            switch ascent                                                                                                   %
                case {'RCGPR'}     %
                    % reset beta after 20 iterations %
                    if mod(ITER,20)==0                              %
                        beta=0;                                     %
                    else                                            %
                        beta = max(real(beta_num/beta_denom), 0);   %
                    end
            end
            %%%
    
            save_gradJ(ITER,1:length(gradJ)) = gradJ;               %
            save_direction(ITER,1:length(direction)) = direction;   %
            %%% Gradient optimization iterative step %%%
            switch ascent                                                                                       %
                case {'RCGPR'}                                                                        %
                    proj_gradJ = projection(gradJ,psi,K0,N,lambda);                                             %
                    psi_x = (2*pi*1i*K0).*psi;                                                                  % derivative of projector term
                    trans_dir = projection(direction,psi,K0,N,lambda) ...                                       %
                        / ( sum(abs(psi).^2) + lambda^2*sum(abs(psi_x).^2) );                                   % vector transport of previous direction
                    direction = proj_gradJ + beta*trans_dir;                                                    % direction towards optimal solution
                    [tau,Jeval,iterJ] = optimize_tau(phi,direction,tau,E,K1,K0,T,nu,N,stepscale,Jeval,iterJ);   % compute optimal step size 
                    phi = phi + tau*direction;                                                                  % iterative solution update
                    phi = retraction(phi,E,K0,N);                                                               % retraction to constraint manifold
                    psi = phi;                                                                                  % projection term
                    fprintf('%01d, %02dh%02dm%02ds\t\t',iterJ,floor(toc/3600),floor(mod(toc/60,60)),floor(mod(toc,60)))
            end
            
            uField_old = uField;                                                %
            Ntime_old = Ntime;                                                  %
            [tvector,uField] = BurgersDS_Fourier(phi,K1,K0,T,nu,N,stepscale);   %
            Ntime = length(tvector);                                            %
            uu = uField(Ntime,:);                                               %
            J_new = eval_J(uu, phi,K0,N);                                       %
            iterJ = iterJ + 1;                                                  %
            Jeval(iterJ) = J_new;                                               %
            fprintf('%01dh%02dm%02ds\t',floor(toc/3600),floor(mod(toc/60,60)),floor(mod(toc,60)))
            
            J_step = (J_new - J(ITER))/J(ITER);                             %
            ITER = ITER+1;                                                  %
            J(ITER) = J_new;                                                %
            K(ITER) = norm(uu);                                             %
            switch ascent                                                                                                   %
                case {'RCGPR'}     %
                    step(ITER) = tau;                                                                                       %
            end
            momentum(ITER) = beta;                                  %
            fprintf('%.4f\t%d\n%02d\t%.4f\t', J(ITER), abs(J_step), ITER, J(ITER))
        end
    
        if abs(J_step) < 0 && ITER > 2      %
            phi = phi_old;                  % cancel out negative step
            uField = uField_old;            %
            Ntime = Ntime_old;              %
        end
    
        %runtime = toc;                                                          % stop timer
        
        %fprintf('P = %d solved' , p );          %
    
        J = rmmissing(J);
    
        disp('Solved optimization problem!')
        fprintf('Iterations \t Initial J \t Optimal J \t Wall Clock \n ')
        fprintf('----------------------------------------------------------------------------------------------------------\n')
        fprintf('%02d \t %.4f\t %f \t ', length(J), J(1), J(end) )
        disp(datetime)
        
        if ITER > MAXITER                                                           %
            disp('Optimization reached maximum number of iterations.');             %
        end
    
        switch VERBOSE                                                                                                  %
            case 1                                                                                                      % save data
                diagnostics_save(J,K,step,momentum,timept,lambda,E0,ITER,testcase);  
                fprintf('Saved diagnostics file at %01dh%02dm%02ds\t \n',floor(toc/3600),floor(mod(toc/60,60)),floor(mod(toc,60))')%
                evaluation_save(Jeval,timept,lambda,E0,iterJ,testcase);                                                 %
                fprintf('Saved objective evaluation file at %01dh%02dm%02ds\t \n',floor(toc/3600),floor(mod(toc/60,60)),floor(mod(toc,60))')%
                f_ens = time_evolution_and_optIC_save(Ntime,phi,timept,E0,lambda,T,N,uField,testcase,K0,stepscale);     %
                fprintf('Saved time evolution file at %01dh%02dm%02ds\t \n',floor(toc/3600),floor(mod(toc/60,60)),floor(mod(toc,60))')%
                branch_save(timept,f_ens,E0,lambda,T,testcase);                                                         %
                fprintf('Saved solution branch file at %01dh%02dm%02ds\t \n',floor(toc/3600),floor(mod(toc/60,60)),floor(mod(toc,60))')%
        end
    else
        disp('Solved forward problem!')
        switch VERBOSE                                                                                                  %
            case 1 
                f_ens = time_evolution_and_optIC_save(Ntime,phi,timept,E0,lambda,T,N,uField,testcase,K0,stepscale);     %
                fprintf('Saved time evolution file at %01dh%02dm%02ds\t \n',floor(toc/3600),floor(mod(toc/60,60)),floor(mod(toc,60))')%
        end
    end
return