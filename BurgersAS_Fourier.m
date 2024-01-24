function solution = BurgersAS_Fourier(uu,K1,K0,T,nu,N,uField,tvector)

%{
K1: fourier space domain for derivative 1 
K0: % fourier space domain
T: length of time window
nu: viscosity coefficient
N: number of grid points
uField: solution to unknown variable
tvector: solution timepoints
%}

    u0_hat = (2*pi*K1).^2.*uu; 
    
    UMAX = max(abs(real(ifft(uField(1,:)))));
    DT = 1/(N*UMAX);
    dt = T/50;
    while dt>DT
        dt = 0.95*dt;
    end
    dt = dt / 50; % adjust time step manually
    N_time = ceil(T/dt);
    tvec = linspace(T,0,N_time);
    tvec = tvec.';
    dt = tvec(1) - tvec(2);
    
    for n=1:N_time-1
        t = tvec(n); 
        uu = interp1(tvector,uField,t);
    
        %==========================================================================
        % Solution to adjoint equation using mixed RK3 and Crank-Nicholson
        vx = (2*pi*1i*K0).*u0_hat; 
        R1 = multiply(uu,vx,'fourier');
        y1_hat = ( (1-nu*(4*dt/15)*(2*pi*K1).^2).*u0_hat + (8*dt/15)*R1 )./( 1+nu*(4*dt/15)*(2*pi*K1).^2 );
        
        vx = (2*pi*1i*K0).*y1_hat;
        R0 = R1;
        R1 = multiply(uu,vx,'fourier'); 
        y2_hat = ( (1-nu*(dt/15)*(2*pi*K1).^2).*y1_hat + (5*dt/12)*R1 - (17*dt/60)*R0 )./( 1+nu*(dt/15)*(2*pi*K1).^2 ); 
        
        vx = (2*pi*1i*K0).*y2_hat;
        R0 = R1;
        R1 = multiply(uu,vx,'fourier'); 
        u1_hat = ( (1-nu*(dt/6)*(2*pi*K1).^2).*y2_hat + (3*dt/4)*R1 - (5*dt/12)*R0 )./( 1+nu*(dt/6)*(2*pi*K1).^2 );
    
        u0_hat=u1_hat;
    end
    
    solution = u1_hat;

return