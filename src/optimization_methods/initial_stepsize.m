function stepsize = initial_stepsize(phi,gradJ,K0,lambda)

    gradJ_x = (2*pi*1i*K0).*gradJ;                                          % derivative of gradient of objective functional
    phi_x = 2*pi*1i*K0.*phi;                                                % derivative of physical data
    stepsize = abs(-0.4* ...
                    real( sum(phi.*conj(gradJ)) + lambda^2*sum(phi_x.*conj(gradJ_x)) ) ... 
                    / ( sum(abs(gradJ).^2) + lambda^2*sum(abs(gradJ_x).^2) ) ...
                    );

return