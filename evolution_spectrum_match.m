for lamtest = 4:5
    for tptest = 19:19
    
        s = 0;
        lam = [ 1, 0.5, 0.1, 1.5, 2];
        lambda = lam(lamtest);
        E0 = 1000;
        N = 2048;
        stepscale = 1;
        if s == 0
            testcase = [ '_' num2str(N) '_LuMax_t' num2str(stepscale) '_fine50' ];
            % testcase = [ '_' num2str(N) '_L1Max_t' num2str(stepscale) '_tp(' num2str(tptest) ')' ];
        else
            testcase = [ '_' num2str(N) '_slingshot' num2str(s) 'Max_t' num2str(stepscale) '_short120' ];
        end
    
        tstart = 1;
        tend = 50;
        
        for timept = tstart:tend
        
            terminal_E0_file = [pwd '/data/time_evolution/terminal' testcase...
                '_E0(' num2str(E0) ')_lambda(' num2str(lambda) ').dat'];
        
            optIC_phys_file = [pwd '/data/spectrum/spectrum_E0_1000/optICphys'...
                testcase '_E0(' num2str(E0) ')_Timept_'...
                num2str(timept) '_lambda(' num2str(lambda) ').dat'];
            optIC_four_file = [pwd '/data/spectrum/spectrum_E0_1000/optICfour'...
                testcase '_E0(' num2str(E0) ')_Timept_'...
                num2str(timept) '_lambda(' num2str(lambda) ').dat'];
            optIC_fourfull_file = [pwd '/data/spectrum/spectrum_E0_1000/optICfourfull'...
                testcase '_E0(' num2str(E0) ')_Timept_'...
                num2str(timept) '_lambda(' num2str(lambda) ').dat'];
            deriv_file = [pwd '/data/spectrum/spectrum_E0_1000/deriv_optICphys'...
                testcase '_E0(' num2str(E0) ')_Timept_'...
                num2str(timept) '_lambda(' num2str(lambda) ').dat'];
            derivfdm_file = [pwd '/data/spectrum/spectrum_E0_1000/derivfdm_optICphys'...
                testcase '_E0(' num2str(E0) ')_Timept_'...
                num2str(timept) '_lambda(' num2str(lambda) ').dat'];
        
            t_evolution = readmatrix(terminal_E0_file);
            optPhi_phys = readmatrix(optIC_phys_file);
            optPhi_four = readmatrix(optIC_four_file);
            optPhi_fourfull = readmatrix(optIC_fourfull_file);
            du_phys = readmatrix(deriv_file);
            du_phys_fdm = readmatrix(derivfdm_file);
        
            tindex = find( isnan(t_evolution(:,2*timept)) , 1, 'first');
            skip = size(optPhi_phys,2) - tindex + 2; 
            if timept == tend
                skip = size(optPhi_phys,2) - size(t_evolution,1) + 1;
            end
    
            if skip > 0 % match up indices by removing insignificant columns
                optPhi_phys = [ optPhi_phys(:,1:10) , optPhi_phys(:,10+skip:end) ]; 
                optPhi_four = [ optPhi_four(:,1:10) , optPhi_four(:,10+skip:end) ]; 
                optPhi_fourfull = [ optPhi_fourfull(:,1:10) , optPhi_fourfull(:,10+skip:end) ]; 
                du_phys = [ du_phys(:,1:10) , du_phys(:,10+skip:end) ]; 
                du_phys_fdm = [ du_phys_fdm(:,1:10) , du_phys_fdm(:,10+skip:end) ]; 
            
                % Save spectrum files
                writematrix(optPhi_phys, optIC_phys_file,'Delimiter','tab');
                writematrix(optPhi_four, optIC_four_file,'Delimiter','tab');
                writematrix(optPhi_fourfull, optIC_fourfull_file,'Delimiter','tab');
                writematrix(du_phys, deriv_file,'Delimiter','tab');
                writematrix(du_phys_fdm, derivfdm_file,'Delimiter','tab');
            end
    
        end
    end
end