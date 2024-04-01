function solution_derivative(E0,filename,TW)

%%% find file %%%
optIC_phys_file = [pwd '/data/spectrum/spectrum_E0_' num2str(E0) '/' filename '.dat'];
optIC_phys = readmatrix(optIC_phys_file);
t = linspace(0,TW,size(optIC_phys,2)-1);
N = size(optIC_phys,1);
K0 = [0:N/2-1 0 -N/2+1:-1];
% x = linspace(1,N,N);
% [T,X] = meshgrid(t,x);
Ntime = length(t);
dt = t(2)-t(1);

uf_x = NaN(N,Ntime+1); % initialize
du = NaN(N,Ntime+1); % initialize
du(:,1) = optIC_phys(:,1); % grid points
uf_x(:,1) = optIC_phys(:,1); % initialize
% differentiate
for i = 2 : Ntime+1

    uf = fft(optIC_phys(:,i));
    uf_x(:,i) = 2*pi*1i*K0'.*uf;    

    du(1,i) = abs( ( optIC_phys(2,i) - 2*optIC_phys(1,i) + optIC_phys(N,i) ) / (dt^2) ) ;
    for j = 2 : N-1
        du(j,i) = abs( ( optIC_phys(j+1,i) - 2*optIC_phys(j,i) + optIC_phys(j-1,i) ) / (dt^2) ) ;
    end
    du(N,i) = abs( ( optIC_phys(1,i) - 2*optIC_phys(N,i) + optIC_phys(N-1,i) ) / (dt^2) ) ;

end

figure(1)
plot(du(:,ceil(Ntime/2)))
figure(2)
plot(abs(real(ifft(uf_x(:,ceil(Ntime/2))))))

% save file
deriv_file = [pwd '/data/spectrum/spectrum_E0_' num2str(E0) '/deriv_' filename '.dat'];
writematrix(du, deriv_file,'Delimiter','tab');
