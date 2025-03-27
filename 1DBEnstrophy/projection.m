function result = projection(original,proj_term,K0,N,lambda)

    original_x = (2*pi*1i*K0).*original; % compute derivative
    proj_term_x = (2*pi*1i*K0).*proj_term; % compute derivative
    numer_proj = sum(proj_term.*conj(original))/N^2 ...
            + lambda^2*sum(proj_term_x.*conj(original_x))/N^2; % < proj, original >
    denomin_proj = sum(abs(proj_term).^2)/N^2 + lambda^2*sum(abs(proj_term_x).^2)/N^2; % < proj, proj >
    projector = real(numer_proj/denomin_proj);
    result = original - projector*proj_term; % Projection of original term

return