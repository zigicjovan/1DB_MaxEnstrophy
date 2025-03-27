function result = multiply(u,v,space)

% Multiply two vector u,v using DEALIASING
% 'space' refers to the space on which the vectors are defined:
% space = 'fourier' or space = 'real'
    
    N = length(u);
    
    n2 = floor(N/2);
    if mod(n2,2)~=0
        n2=n2-1;
    end
    
    myfilter = [ones(1,n2/2) zeros(1,N-n2) ones(1,n2/2)]; 
    
    switch space
        case 'real'
            u_hat = myfilter.*fft(u);
            v_hat = myfilter.*fft(v);
    
            u = real(ifft(u_hat));
            v = real(ifft(v_hat));
            
            aux = myfilter.*fft(u.*v);
            
            result = real(ifft(aux));
            
        case 'fourier'
            u_hat = myfilter.*u;
            v_hat = myfilter.*v;
            
            u = real(ifft(u_hat));
            v = real(ifft(v_hat));
            result = myfilter.*fft(u.*v);
    end

return