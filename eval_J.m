function solution = eval_J(u,phi,K0,N)

    u_x = 2*pi*1i*K0.*u;
    phi_x = 2*pi*1i*K0.*phi;
    solution = 0.5*( sum(abs(u_x).^2) - sum(abs(phi_x).^2) )/N^2; 

return