function result = optimize_tau(phi,grad,iniTau,Ens,K1,K0,T,nu,N)
% find optimal stepsize for objective functional using Brent's method %

    VERBOSE = 0;
    
    TOL = 10^(-5);
    GOLD = 1.618034; GLIMIT = 100;
    ITMAX = 300; CGOLD = .381966; ZEPS = 10^(-10);
    
    % initialize 3 x-points (left, center, right) 
    % and keep switching them to converge on minimum value
    
    % This part of the code brackets the location of the minimum of the
    % functional J
    iter=0;
    tA = 0; % initial left x
    tB = iniTau + 10^(-9); % initial center x
    
    phi_bar = phi + tA*grad;
    phi_bar = adjust_optIC(phi_bar,Ens,K0,N);
    [tvec,u] = BurgersDS_Fourier(phi_bar,K1,K0,T,nu,N);
    Ntime = length(tvec);
    uu = u(Ntime,:);
    FA = -eval_J(uu,phi_bar,K0,N); % initial left f(x)
    
    if VERBOSE == 5
        tau_prog(tA, FA); % tau and f(tau) % checkpoint 0
    end
        %brak_prog(uu, tA, tB, phi, phi_bar, grad, Ens); % checkpoint 0
    
    phi_bar = phi + tB*grad;
    phi_bar = adjust_optIC(phi_bar,Ens,K0,N);
    [tvec,u] = BurgersDS_Fourier(phi_bar,K1,K0,T,nu,N);
    Ntime = length(tvec);
    uu = u(Ntime,:);
    FB = -eval_J(uu,phi_bar,K0,N); % initial center f(x)
    
    % tau_prog(tB, FB); % tau and f(tau) % checkpoint 1
    if VERBOSE == 5
        tau_prog(tB, FB); 
    end
    
    if FB>FA %% existing code leads to negative tC value
        aux = tA;
        tA = tB; % new center x
        tB = aux; % new right x
        aux = FA;
        FA = FB; % new center f(x)
        FB = aux; % new right f(x)
    end
    
    % FA > FB, tB ? tA
    
    tC = tB + GOLD*(tB - tA); % initial right or center x %% abs to stay positive??
    phi_bar = phi + tC*grad;
    phi_bar = adjust_optIC(phi_bar,Ens,K0,N);
    [tvec,u] = BurgersDS_Fourier(phi_bar,K1,K0,T,nu,N);
    Ntime = length(tvec);
    uu = u(Ntime,:);
    FC = -eval_J(uu,phi_bar,K0,N); % initial right or center f(x)
    
    % brak_prog(uu, tA, tC, phi, phi_bar, grad, Ens); % checkpoint 1
    % brent_prog(tA, tB, FB, tC, phi, phi_bar, grad, Ens); % checkpoint 1
    
    while FB>=FC && iter<ITMAX 
        
        %     SA = (tC-tB)*FA;
        %     SB = (tA-tC)*FB;
        %     SC = (tB-tA)*FC;
        %     tP = 0.5*( (tC+tB)*SA + (tA+tC)*SB + (tB+tA)*SC)/(SA+SB+SC);
        iter = iter+1;
        R = (tB-tA)*(FB-FC); 
        Q = (tB-tC)*(FB-FA); 
        tP = tB - 0.5*((tB-tC)*Q - (tB-tA)*R)/(sign(Q-R)*abs(Q-R)); % new x correction
        
        Pmax = tB + GLIMIT*(tC-tB); %% abs to stay positive??
    
        if (tB-tP)*(tP-tC)>0 % if tP is center x
            phi_bar = phi + tP*grad;
            phi_bar = adjust_optIC(phi_bar,Ens,K0,N);
            [tvec,u] = BurgersDS_Fourier(phi_bar,K1,K0,T,nu,N);
            Ntime = length(tvec);
            uu = u(Ntime,:);
            FP = -eval_J(uu,phi_bar,K0,N); % new f(x)
    
            % tau_prog(tP, FP); % tau and f(tau) % checkpoint 2
            if VERBOSE == 5
                tau_prog(tP, FP); 
            end
            %brak_prog(uu, tB, tC, phi, phi_bar, grad, Ens); % checkpoint 2
            
            if FP<FC
                tA = tB; % new end x
                FA = FB; % new f(x)
                tB = tP; % new center x
                FB = FP; % new f(x)
                break;
            elseif FP>FB
                tC = tP;
                FC = FP;
                break;
            end
            
            tP = tC + GOLD*(tC-tB); % new x correction
            phi_bar = phi + tP*grad;
            phi_bar = adjust_optIC(phi_bar,Ens,K0,N);
            [tvec,u] = BurgersDS_Fourier(phi_bar,K1,K0,T,nu,N);
            Ntime = length(tvec);
            uu = u(Ntime,:);
            FP = -eval_J(uu,phi_bar,K0,N);
    
            % tau_prog(tP, FP); % tau and f(tau) % checkpoint 3
            if VERBOSE == 5
                tau_prog(tP, FP); 
            end
            %brak_prog(uu, tA, tC, phi, phi_bar, grad, Ens); % checkpoint 3
            
        elseif (tC-tP)*(tP-Pmax)>0
            phi_bar = phi + tP*grad;
            phi_bar = adjust_optIC(phi_bar,Ens,K0,N);
            [tvec,u] = BurgersDS_Fourier(phi_bar,K1,K0,T,nu,N);
            Ntime = length(tvec);
            uu = u(Ntime,:);
            FP = -eval_J(uu,phi_bar,K0,N);
            
            if FP<FC
                tB = tC;
                tC = tP;
                FB = FC;
                FC = FP;
                tP = tC+GOLD*(tC-tB); % new x correction
                phi_bar = phi + tP*grad;
                phi_bar = adjust_optIC(phi_bar,Ens,K0,N);
                [tvec,u] = BurgersDS_Fourier(phi_bar,K1,K0,T,nu,N);
                Ntime = length(tvec);
                uu = u(Ntime,:);
                FP = -eval_J(uu,phi_bar,K0,N);
            end
    
            % tau_prog(tP, FP); % tau and f(tau) % checkpoint 4
            if VERBOSE == 5
                tau_prog(tP, FP); 
            end
            %brak_prog(uu, tA, tC, phi, phi_bar, grad, Ens); % checkpoint 4
            
        elseif (tP-Pmax)*(Pmax-tC)>=0
            tP = Pmax; % new x correction
            phi_bar = phi + tP*grad;
            phi_bar = adjust_optIC(phi_bar,Ens,K0,N);
            [tvec,u] = BurgersDS_Fourier(phi_bar,K1,K0,T,nu,N);
            Ntime = length(tvec);
            uu = u(Ntime,:);
            FP = -eval_J(uu,phi_bar,K0,N);
    
            if VERBOSE == 5
                tau_prog(tP, FP); % tau and f(tau) % checkpoint 5
            end
            %brak_prog(uu, tA, tC, phi, phi_bar, grad, Ens); % checkpoint 5
            
        else
            tP = tC + GOLD*(tC-tB); % new x correction
            phi_bar = phi + tP*grad;
            phi_bar = adjust_optIC(phi_bar,Ens,K0,N);
            [tvec,u] = BurgersDS_Fourier(phi_bar,K1,K0,T,nu,N);
            Ntime = length(tvec);
            uu = u(Ntime,:);
            FP = -eval_J(uu,phi_bar,K0,N);
    
            % tau_prog(tP, FP); % tau and f(tau) % checkpoint 6
            if VERBOSE == 5
                tau_prog(tP, FP); 
            end
            %brak_prog(uu, tA, tC, phi, phi_bar, grad, Ens); % checkpoint 6
        end
        
        tA = tB;
        tB = tC;
        tC = tP;
        FA = FB;
        FB = FC;
        FC = FP;
    
    end
    
    if iter==ITMAX
        result = iniTau;
        return
    end
    
    % This part of the code finds the minimum of the functional and gives
    % the value of the optimal tau
    
    % X = run_brent(tA, tB, tC, FA, FB, FC, TOL, ZEPS, CGOLD, ITMAX, phi, grad, Ens);

    D = 0;
    A = min(tA,tC);
    B = max(tA,tC);
    V = tB; W = V; X = V; E = 0;
    FX = FB; FV = FX; FW = FX;
    % U_tol = 10^(-9.71);
    % M = 10;
  
    for j=1:ITMAX
              
        XM = 0.5*(A+B); %halfway between [A,B]
        TOL1 = (TOL*abs(X)+ZEPS);
        TOL2 = 2*TOL1;
        x_step = abs(X-XM);
        x_tol = (TOL2-0.5*(B-A));
    
        if ( x_step <= x_tol )
            break;
        end

        FLAG = 1;
        if ( abs(E) > TOL1 )
            R = (X-W)*(FX-FV);
            Q = (X-V)*(FX-FW);
            P = (X-V)*Q - (X-W)*R;
            Q = 2*(Q-R);
            if Q>0
                P = -P;
            end
            Q = abs(Q);
            ETEMP = E;
            E = D;
            
            if (abs(P) >= abs(0.5*Q*ETEMP)) || (P <= Q*(A-X)) || ...
                    (P >= Q*(B-X))
                FLAG = 1;
            else
                FLAG = 2;
            end
            
        end
        
        switch FLAG
            case 1
                if X >= XM
                    E = A-X;
                else
                    E=B-X;
                end
                D = CGOLD*E;
            case 2
                D = P/Q;
                U = X+D;
                if (U-A < TOL2) || (B-U < TOL2)
                    D = sign(XM-X)*TOL1;
                end
        end
        
        if abs(D) >= TOL1
            U = X+D;
        else
            U = X + sign(D)*TOL1;
        end
        
        phi_bar = phi + U*grad;
        phi_bar = adjust_optIC(phi_bar,Ens,K0,N);
        [tvec,u] = BurgersDS_Fourier(phi_bar,K1,K0,T,nu,N);
        Ntime = length(tvec);
        uu = u(Ntime,:);
        FU = -eval_J(uu,phi_bar,K0,N);

        % tau_prog(U, FU); % tau and f(tau) % checkpoint 7
        if VERBOSE == 5
            tau_prog(U, FU); 
        end
        
        if FU <= FX
            if U >= X
                A = X;
            else
                B = X;
            end
            V = W;
            FV = FW;
            W = X;
            FW = FX;
            X = U;
            FX = FU;
        else
            if U < X
                A = U;
            else
                B = U;
            end
            
            if (FU <= FW) || (W == X)
                V = W;
                FV = FW;
                W = U;
                FW = FU;
            elseif (FU <= FV) || (V==X) || (V==W)
                V = U;
                FV = FU;
            end
        end
    
    end

    result = X;

return