function solution = eval_grad_J(u_star,phi,K1,lambda)

    phi_xx = -(2*pi*K1).^2.*phi;
    aux = u_star + phi_xx; 
    solution = aux./(1 + (lambda*2*pi*K1).^2);

return