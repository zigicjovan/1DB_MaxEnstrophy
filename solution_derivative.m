function solution_derivative(E0,filename,TW)

%%% find file %%%
optIC_phys_file = [pwd '/data/spectrum/spectrum_E0_' num2str(E0) '/' filename '.dat'];
optIC_phys = readmatrix(optIC_phys_file);
t = linspace(0,TW,size(optIC_phys,2)-1);
N = size(optIC_phys,1);
% x = linspace(1,N,N);
% [T,X] = meshgrid(t,x);
Ntime = length(t);
dt = t(2)-t(1);

du = NaN(N,Ntime+1); % initialize
du(:,1) = optIC_phys(:,1); % grid points
% differentiate
for i = 2 : Ntime+1

    % initial term, forwards ( xf - xc / h )
    du(1,i) = abs( ( optIC_phys(2,i) - optIC_phys(1,i) ) / dt ) ;

    % inner terms, centrally ( xf - 2xc + xb / h^2 )
    for j = 2 : N-1
        du(j,i) = abs( ( optIC_phys(j+1,i) - 2*optIC_phys(j,i) + optIC_phys(j-1,i) ) / (dt^2) ) ;
    end
    
    % final term, backwards ( xc - xb / h )
    du(N,i) = abs( ( optIC_phys(N,i) - optIC_phys(N-1,i) ) / dt ) ;

end

% save file
deriv_file = [pwd '/data/spectrum/spectrum_E0_' num2str(E0) '/deriv_' filename '.dat'];
writematrix(du, deriv_file,'Delimiter','tab');
