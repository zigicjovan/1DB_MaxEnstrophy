function [phi,f_ens] = optimize_Burgers(opt_sol,CONT,testcase,enspt,timept,phi,N,T,...
                                        x,Enstrophy,TimeWindow,stepscale,VERBOSE,...
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
    
    fprintf('Solving 1D Burgers, N = %d, E0 = %.0f (point %d of %d), T = %.4f (point %d of %d) \n', ...
        N , E0 , enspt , length(Enstrophy) , T , timept , length(TimeWindow) );        %

    [tvector,uField,Ntime,u_TC] = BurgersDS_Fourier(phi,K1,K0,T,nu,N,stepscale);       % solve forward PDE

    J(ITER) = eval_J(u_TC,phi,K0,N);                                        % evaluate objective functional (E(T) - E0)
    iterJ = iterJ + 1;                                                      %
    Jeval(iterJ) = J(ITER);                                                 %
    K(ITER) = norm(u_TC);                                                   % kinetic energy at E(T)
    step(ITER) = 0;                                                         %
    momentum(ITER) = 0;                                                     %
    save_gradJ = NaN(1,N);
    save_direction = zeros(1,N);                                            %
    J_step = 1;

    if opt_sol == 1

        fprintf('Iter \t J0 \t Adjoint solved \t Line-search \t Forward solved \t J1 \t\t J change \n ')
        fprintf('----------------------------------------------------------------------------------------------------------\n')    
        fprintf('%02d\t%.4f\t', ITER, J(ITER))
    
        while ( abs(J_step) > J_step_tol ) && ( ITER <= MAXITER )           %
            
            u_adj = BurgersAS_Fourier(u_TC,K1,K0,T,nu,N,uField,tvector,stepscale);    % solve adjoint PDE
            gradJ = eval_grad_J(u_adj,phi,K1,lambda);                           %
            if ITER == 1
                tau = initial_stepsize(phi,gradJ,K0,lambda);                          % initial guess for optimal step (tested; do not change)
                proj_gradJ = 1;
                eta = 1;
            end
            beta = momentum_term(ITER,proj_gradJ,gradJ,phi,psi,eta,K0,N,lambda);
            
            fprintf('%01dh%02dm%02ds\t',floor(toc/3600),floor(mod(toc/60,60)),floor(mod(toc,60)))
    
            save_gradJ(ITER,1:length(gradJ)) = gradJ;               %
            
            %% Gradient optimization iterative step                                  %
            proj_gradJ = projection(gradJ,psi,K0,N,lambda);                                             %
            if ITER == 1
                direction = proj_gradJ;                                                                 % direction towards optimal solution
            else
                save_direction(ITER,1:length(direction)) = direction;   %
                trans_dir = vectortransport(phi,direction,eta,K0,N,lambda);                             % vector transport of previous direction
                direction = proj_gradJ + beta*trans_dir;                                                % direction towards optimal solution
            end
            [tau,Jeval,iterJ] = optimize_stepsize(phi,direction,tau,E,K1,K0,T,nu,N,stepscale,Jeval,iterJ);   % compute optimal step size 
            eta = phi + tau*direction;                                                                  % iterative solution update
            phi = retraction(eta,E,K0,N);                                                               % retraction to constraint manifold
            psi = phi;                                                                                  % projection term
            %%
            fprintf('%01d, %02dh%02dm%02ds\t\t',iterJ,floor(toc/3600),floor(mod(toc/60,60)),floor(mod(toc,60)))
            
            [tvector,uField,Ntime,u_TC] = BurgersDS_Fourier(phi,K1,K0,T,nu,N,stepscale);   %
  
            J_new = eval_J(u_TC,phi,K0,N);                                       %
            iterJ = iterJ + 1;                                                  %
            Jeval(iterJ) = J_new;                                               %
            
            fprintf('%01dh%02dm%02ds\t',floor(toc/3600),floor(mod(toc/60,60)),floor(mod(toc,60)))
            
            J_step = (J_new - J(ITER))/J(ITER);                             %
            ITER = ITER + 1;                                                  %
            J(ITER) = J_new;                                                %
            K(ITER) = norm(u_TC);     
            step(ITER) = tau;                                                                                       %
            momentum(ITER) = beta;                                  %
            
            fprintf('%.4f\t%d\n%02d\t%.4f\t', J(ITER), abs(J_step), ITER, J(ITER))

            switch VERBOSE                                                                                          %
                case 3                                                                                              % test definition of Gateaux differential
                    u_adj = BurgersAS_Fourier(u_TC,K1,K0,T,nu,N,uField,tvector,stepscale);                            %
                    gradJ = eval_grad_J(u_adj,phi,K1,lambda);                                                   %
                    gradJ_x = (2*pi*1i*K0).*gradJ;                                                          %
                    [epsilon,kappa] = kappaTestFourier(phi,gradJ,gradJ_x,J(ITER),x,K1,K0,T,nu,N,lambda,stepscale);  %
                    kappa_save(kappa, epsilon, timept,E0,lambda,ITER,testcase);                                     %
            end
        end
    
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