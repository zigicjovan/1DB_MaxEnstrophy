function solution_in_time_figure(E0,filename,TW,casename)

close all

%%%% adjust parameters %%%
% E0 = 1000;
% filename = [ 'optICphys_2048_slingshot5Tfix_long27_E0(' num2str(E0) ')_Timept_14_lambda(0.5)' ];
% TW = 0.028563;

%%% find file %%%
optIC_phys_file = [pwd '/data/spectrum/spectrum_E0_' num2str(E0) '/' filename '.dat'];
optIC_phys = readmatrix(optIC_phys_file);
t = linspace(0,TW,size(optIC_phys,2)-1);
N = size(optIC_phys,1);
x = linspace(1,N,N);
[T,X] = meshgrid(t,x);
Ntime = length(t);

%%% make gif %%%

%
for i = 2 : ceil(Ntime/500) : Ntime
    pbaspect( [ max(max(t)), max(max(x)), max(max(optIC_phys(:,2:end))) ] );
    surf(T(:,1:i),X(:,1:i),optIC_phys(:,2:(2+i-1)))
    shading(gca,'interp')
    colormap(redblue)
    xlabel('Time')
    ylabel('Space')
    zlabel('Solution')
    title(['1DB, E0 = ' num2str(E0) ', T = ' num2str(TW) ', N = ' num2str(N) ''])
    axis( [ 0 max(max(t))  0 max(max(x))  min(min(optIC_phys(:,2:(2+i-1)))) max(max(optIC_phys(:,2:(2+i-1)))) ] )
    view([37.5,30]);
    drawnow

    if i == 2
        gif([ 'solution1DB_' casename '_evolution_E0' num2str(E0) '_N' num2str(N) '_T' num2str(TW) '_3Dview.gif'])
    else
        gif
    end
end
close all
%}
for i = 2 : ceil(Ntime/500) : Ntime
    pbaspect( [ max(max(t)), max(max(x)), max(max(optIC_phys(:,2:end))) ] );
    surf(T(:,1:i),X(:,1:i),optIC_phys(:,2:(2+i-1)))
    shading(gca,'interp')
    colormap(redblue)
    xlabel('Time')
    ylabel('Space')
    zlabel('Solution')
    title(['1DB, E0 = ' num2str(E0) ', T = ' num2str(TW) ', N = ' num2str(N) ''])
    axis( [ 0 max(max(t))  0 max(max(x))  min(min(optIC_phys(:,2:(2+i-1)))) max(max(optIC_phys(:,2:(2+i-1)))) ] )
    view([0,-45,0]);
    drawnow

    if i == 2
        gif([ 'solution1DB_' casename '_evolution_E0' num2str(E0) '_N' num2str(N) '_T' num2str(TW) '_sideview.gif'])
    else
        gif
    end
end
close all
%}
for i = 2 : ceil(Ntime/500) : Ntime
    pbaspect( [ max(max(t)), max(max(x)), max(max(optIC_phys(:,2:end))) ] );
    surf(T(:,1:i),X(:,1:i),optIC_phys(:,2:(2+i-1)))
    shading(gca,'interp')
    colormap(redblue)
    xlabel('Time')
    ylabel('Space')
    zlabel('Solution')
    title(['1DB, E0 = ' num2str(E0) ', T = ' num2str(TW) ', N = ' num2str(N) ''])
    axis( [ 0 max(max(t))  0 max(max(x))  min(min(optIC_phys(:,2:(2+i-1)))) max(max(optIC_phys(:,2:(2+i-1)))) ] )
    view(2);
    drawnow

    if i == 2
        gif([ 'solution1DB_' casename '_evolution_E0' num2str(E0) '_N' num2str(N) '_T' num2str(TW) '_topview.gif'])
    else
        gif
    end
end
%}

