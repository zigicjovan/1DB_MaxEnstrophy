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

    % Ntime = length(tvector);
    % uu = uField(Ntime,:);
    u0_hat = (2*pi*K1).^2.*uu; 
    
    % factor = 0.00025;%tmp
    % filter = 1./(1+(factor*2*pi*K1).^2);%tmp
    
    UMAX = max(abs(real(ifft(uField(1,:)))));
    DT = 1/(N*UMAX);
    dt = T/50;
    while dt>DT
        dt = 0.95*dt;
    end
    % dt = dt / 10; %[advecfix5: dt/2], [advecfix6: dt/10]
    N_time = ceil(T/dt);
    tvec = linspace(T,0,N_time);
    tvec = tvec.';
    dt = tvec(1) - tvec(2);
    
    for n=1:N_time-1
        t = tvec(n); %changed from n to N_time - n [advecfix1]
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
    
        % u1_hat = filter.*u1_hat; %tmp [advecfix2] [comment out - advecfix4]
    
        u0_hat=u1_hat;
    end
    
    solution = u1_hat;

return

%{
if toPlot
    plot(x,real(ifft(solution)));
    legend('Initial Condition', 'Solution at t = 0');
    axis([0 1 -5000 5000]);
end

function sol = RHS(t,y)

global nu; 
global uField; global tvector;
global K0; global K1; global I;

% I = sqrt(-1);
% K0 = [0:N/2-1 0 -N/2+1:-1];
% K1 = [0:N/2 -N/2+1:-1];

vx = (2*pi*I*K0).*y.';
uu = interp1(tvector,uField,t,'spline');
sol = nu*(2*pi*K1).^2.*y.' - multiply(uu,vx,'fourier');
sol = sol.';

return


function sol = RHSjacobian(t,y)

global nu; global N; 
global uField; global tvector;
global K0; global K1; global I;

sol = zeros(N,N);
uu = interp1(tvector,uField,t,'spline');
uu = fliplr(uu);

for j=1:N
    sol(j,:) = uu.*(2*pi*I*K0);
    uu = circshift(uu, [0 1]);
end

sol = sol - nu*diag((2*pi*K1).^2);
return
%}