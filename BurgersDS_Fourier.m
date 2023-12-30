function [tvec,solution]=BurgersDS_Fourier(u0_hat,K1,K0,T,nu,N)

%{
K1: fourier space domain for derivative 1 
K0: fourier space domain
T: length of time window
nu: viscosity coefficient
N: number of grid points
%}
    
    % factor = 0.00025;
    % filter = 1./(1+(factor*2*pi*K1).^2);
    
    UMAX = max(abs(real(ifft(u0_hat))));
    DT = 1/(N*UMAX); %timestep below this
    dt = T/50;
    while dt>DT
        dt = 0.95*dt;
    end
    % dt = dt / 10; %[advecfix5: dt/2], [advecfix6: dt/10]
    Ntime = ceil(T/dt);
    tvec = linspace(0,T,Ntime);
    tvec = tvec.';
    dt = tvec(2) - tvec(1);
    
    solution = zeros(Ntime,N);
    solution(1,:) = u0_hat;
    
    for i=2:Ntime
        
        % Solution to Burgers equation using mixed RK3 and Crank-Nicholson
        R1 = -0.5*(2*pi*1i*K0).*multiply(u0_hat,u0_hat,'fourier'); 
        y1_hat = ( (1-nu*(4*dt/15)*(2*pi*K1).^2).*u0_hat + (8*dt/15)*R1 )./( 1+nu*(4*dt/15)*(2*pi*K1).^2 );
        
        R0 = R1;
        R1 = -0.5*(2*pi*1i*K0).*multiply(y1_hat,y1_hat,'fourier');
        y2_hat = ( (1-nu*(dt/15)*(2*pi*K1).^2).*y1_hat + (5*dt/12)*R1 - (17*dt/60)*R0 )./( 1+nu*(dt/15)*(2*pi*K1).^2 ); 
        
        R0 = R1;
        R1 = -0.5*(2*pi*1i*K0).*multiply(y2_hat,y2_hat,'fourier');
        u1_hat = ( (1-nu*(dt/6)*(2*pi*K1).^2).*y2_hat + (3*dt/4)*R1 - (5*dt/12)*R0 )./( 1+nu*(dt/6)*(2*pi*K1).^2 );
    
        % u1_hat = filter.*u1_hat; % [advecfix3 - no filter]
        solution(i,:) = u1_hat;
        u0_hat = u1_hat;
            
    end

return

%{
if toPlot
    figure; hold on; grid on;
    
    plot(x,real(ifft(solution(1,:))),'--r');
    plot(x,real(ifft(solution(Ntime,:))));
    xlabel('x')
    ylabel('u(x)');
    nuVal = num2str(nu);
    myTitle = ['Solution to Burgers Equation with \nu = ' nuVal];
    title(myTitle);
    legend('Initial condition','Solution at time T');
end

function sol = RHS(t,y)

global nu;
global K0; global K1; global I;

sol = -nu*(2*pi*K1).^2.*y.' - 0.5*(2*pi*I*K0).*multiply(y.',y.','fourier');
sol = sol.';

return

function sol = RHSjacobian(t,y)

global nu; global K1;

vec = -nu*(2*pi*K1).^2;
sol = diag(vec);

return
%}