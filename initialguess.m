function solution = initialguess(type,x,E0,a1,x_ref,timept,time_start)
    
    switch type
        case 'sine'
            solution = sqrt(E0)/pi*sin(2*pi*x - pi);
        case 'sin2'
            solution = sqrt(E0)/(2*pi)*sin(4*pi*x - pi);
        case 'exact'
            initialdata_file = [pwd '/optIC/LuLu_u_E010_nu3_N2048.dat'];
            aux = readmatrix(initialdata_file);
            if length(x) == 2048
                solution = (aux).';
            else
                solution = interp1(x_ref, aux , x, 'spline');
            end
            solution = fft(solution);
        case 'optIC'
            try
                phi_guess_file = [pwd '/optIC/phi_E0_' num2str(E0) '__4096_optIC2048ContLu7.dat'];
                phi_guess = readmatrix(phi_guess_file);
                solution = phi_guess(timept,:);
                solution = real(ifft(solution));
                length_ref = length(solution);
                if length(x) ~= length_ref
                    solution = interp1(x_ref, solution , x, 'spline');
                end
                solution = fft(solution);
            catch
                disp('No optIC file found, using exact solution');
                phi = initialguess('exact', x, E0, 0, x_ref);
                solution = fft(phi);
            end
        case '2 modes'
            solution = a1*sin(2*pi*x-pi) + sqrt( E0-(a1*pi)^2 )/(2*pi)*sin(4*pi*x-pi/7);
    end
return