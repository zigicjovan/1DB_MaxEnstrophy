function Du = derivative(u,order)

    numPoints = length(u);
    K1 = [0:numPoints/2-1 0 -numPoints/2+1:-1];
    K2 = [0:numPoints/2 -numPoints/2+1:-1];
    I = sqrt(-1);
    
    parity = mod(order,4);
    switch parity
        case 1
            Du_hat = I*(2*pi*K1).^(order).*fft(u);
            Du = ifft(Du_hat);
        case 2
            Du_hat = -(2*pi*K2).^(order).*fft(u);
            Du = ifft(Du_hat);
        case 3
            Du_hat = -I*(2*pi*K1).^(order).*fft(u);
            Du = ifft(Du_hat);
        case 0
            Du_hat = (2*pi*K2).^(order).*fft(u);
            Du = ifft(Du_hat);
    end
    
    % figure; hold on;
    % plot(u,'r');
    % plot(Du,'.');

return

