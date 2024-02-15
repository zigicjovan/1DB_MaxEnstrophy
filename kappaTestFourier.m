function [epsilon,kappa] = kappaTestFourier(phi,grad,grad_x,J_phi,x,K1,K0,T,nu,N,lambda,stepscale)

%{
phi: optimal IC
grad: gradient of objective functional
grad_x: derivative of grad
J_phi: objective functional
x: physical space domain
K0: fourier space domain
N: number of grid points
lambda: length-scale parameter
%}

    NumPoints = 30;
    NumModes = 1;
    
    epsilon = logspace(-15,0,NumPoints);
    kappa = zeros(length(epsilon),NumModes);
    numerator = zeros(length(epsilon),NumModes);
    denominator = zeros(length(epsilon),NumModes);
    modes = 1:NumModes;

    for k=1:NumModes
        
        phi_pert = 0.5*sin(2*pi*modes(k)*(x - 0.05));
        phi_pert = fft(phi_pert);
        phi_pert_x = (2*pi*1i*K0).*phi_pert;
        innerProd = real(sum(grad.*conj(phi_pert))/N^2 + lambda^2*sum(grad_x.*conj(phi_pert_x))/N^2);
    
        for j=1:NumPoints
            phiBar = phi + epsilon(j)*phi_pert;
            
            [tvec,u] = BurgersDS_Fourier(phiBar,K1,K0,T,nu,N,stepscale);
            Ntime = length(tvec);
            uu = u(Ntime,:);
            ux = (2*pi*1i*K0).*uu;
            phiBar_x = (2*pi*1i*K0).*phiBar;
            J_phiBar = 0.5*sum(abs(ux).^2)/N^2 - 0.5*sum(abs(phiBar_x).^2)/N^2; 
    
            numerator(j,k) = (J_phiBar - J_phi)/(epsilon(j));
            denominator(j,k) = innerProd;
    
            kappa(j,k) = (J_phiBar - J_phi)/(epsilon(j)*innerProd);
        end
        
        % kappa_out = [ epsilon', numerator(:,k), denominator(:,k), kappa(:,k)]
        
    end

return