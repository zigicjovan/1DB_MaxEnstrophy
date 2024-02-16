close all

optIC_phys_file = [pwd '/data/spectrum/spectrum_E0_31622.7766/optICphys_6144_optIC2048ContLu7_t5_E0(31622.7766)_Timept_12_lambda(0.5)_case3.dat'];
optIC_phys = readmatrix(optIC_phys_file);
t = linspace(0,0.00435361025953819,1717);
x = linspace(1,6144,6144);
[T,X] = meshgrid(t,x);
%{
surf(T(:,1:end),X(:,1:end),optIC_phys(:,9:end-1))
shading(gca,'interp')
axis( [ 0 max(t)  0 max(x)  min(optIC_phys(:)) -min(optIC_phys(:)) ] )
xlabel('Time')
ylabel('Space')
zlabel('Solution')
title('1D Burgers, Initial Enstrophy = 31623, Time Window = 0.0043536, Grid points = 6144')
%}

h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
mkdir([pwd  '/data/media/movies' ]);
% filename = [pwd '/data/media/movies/phys_' method '_N_' num2str(N) ...
    % '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1) '_Ls2_' num2str(L_s2) '.gif'];

Ntime = length(t);

% figure(2);
xlabel('Time')
ylabel('Space')
zlabel('Solution')
title('1D Burgers, Initial Enstrophy = 31623, Time Window = 0.0043536, Grid points = 6144')
for i = 2 : ceil(Ntime/500) : Ntime
    % Draw surface plot
    pbaspect( [ max(max(t)), max(max(x)), max(max(optIC_phys(:,9:end-1))) ] );
    surf(T(:,1:i),X(:,1:i),optIC_phys(:,9:(9+i-1)))
    shading(gca,'interp')
    axis( [ 0 max(max(t))  0 max(max(x))  min(min(optIC_phys(:,9:(9+i-1)))) max(max(optIC_phys(:,9:(9+i-1)))) ] )
    %{
    if i < 280
        axis( [ 0 max(max(t))  0 max(max(x))  min(min(optIC_phys(:,9:(9+280-1)))) max(max(optIC_phys(:,9:(9+280-1)))) ] )
    elseif i < 580   
        axis( [ 0 max(max(t))  0 max(max(x))  min(min(optIC_phys(:,9:(9+580-1)))) max(max(optIC_phys(:,9:(9+580-1)))) ] )
    elseif i < 1080
        axis( [ 0 max(max(t))  0 max(max(x))  min(min(optIC_phys(:,9:(9+1080-1)))) max(max(optIC_phys(:,9:(9+1080-1)))) ] )
    elseif i < 1480
        axis( [ 0 max(max(t))  0 max(max(x))  min(min(optIC_phys(:,9:(9+1480-1)))) max(max(optIC_phys(:,9:(9+1480-1)))) ] )
    else 
        axis( [ 0 max(max(t))  0 max(max(x))  min(min(optIC_phys(:,9:end-1))) max(max(optIC_phys(:,9:end-1))) ] )
    end
    %}
    view(3);
    drawnow

    if i == 2
        gif('myfile2.gif')
    else
        gif
    end
    %{
    % Save image
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);

    % Write to GIF File
    if i == 1
        imwrite(imind,cm,'solution1DB','gif', 'DelayTime',0.02, 'Loopcount',inf);
    else
        imwrite(imind,cm,'solution1DB','gif','WriteMode','append');
    end
    %}
end