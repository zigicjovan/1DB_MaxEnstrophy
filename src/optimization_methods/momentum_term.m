function beta = momentum_term(ITER,proj_gradJ,gradJ,phi,psi,eta,K0,N,lambda)

    proj_gradJ_new = projection(gradJ,psi,K0,N,lambda);                 %
    proj_gradJ_x = (2*pi*1i*K0).*proj_gradJ;                                          % derivative of gradient of objective functional
    beta_denom = sum(abs(proj_gradJ).^2)/N^2 + lambda^2*sum(abs(proj_gradJ_x).^2)/N^2;

    proj_gradJ_x_new = (2*pi*1i*K0).*proj_gradJ_new;                        %
    trans_proj_gradJ = vectortransport(phi,proj_gradJ,eta,K0,N,lambda);     % vector transport of previous direction
    trans_proj_gradJ_x = (2*pi*1i*K0).*trans_proj_gradJ;                    %
    beta_num = sum(proj_gradJ_new.*conj(proj_gradJ_new - trans_proj_gradJ))/N^2 ...   %
        + lambda^2*sum(proj_gradJ_x_new.*conj(proj_gradJ_x_new - trans_proj_gradJ_x))/N^2;
               
    % reset beta after 20 iterations %
    if ( mod(ITER,20) == 0 || ITER == 1 )           %
        beta = 0;                                     %
    else                                            %
        beta = max(real(beta_num/beta_denom), 0);   %
    end

return