function [testcase_max,testcase1,testcase0] = testcase_name(caseID,N,stepscale,ascent,tptest,pm,s,shiftT,ss_shift,sstest,Lip,cond,s_w,testnumber)

    if s == 0
        switch ascent
            case {'RHB','RAHB','RNAG','RSNAG','RBNAG','RBNAGa'}
                testcase_max = [ '_' num2str(N) '_LuMax_t' num2str(stepscale) '_' ascent '_L' num2str(Lip) '_K' num2str(cond) '' caseID ];%'p' num2str(caseID) '' ];
                testcase1 = [ '_' num2str(N) '_LuCont1_t' num2str(stepscale) '_' ascent '_L' num2str(Lip) '_K' num2str(cond) '' caseID ];%'p' num2str(caseID) '' ];
                testcase0 = [ '_' num2str(N) '_LuCont0_t' num2str(stepscale) '_' ascent '_L' num2str(Lip) '_K' num2str(cond) '' caseID ];%'p' num2str(caseID) '' ];
            case {'RCGPR','RCGRMIL','RSCG','RMSCG','RDSCG','RJSCG','RSWMN','RSWM','RSWCGN','RSWCGJ','RSWGN','RSWG'}
                testcase_max = [ '_' num2str(N) '_LuMax_t' num2str(stepscale) '_' ascent '' s_w '' caseID '_' testnumber ];
                testcase1 = [ '_' num2str(N) '_LuCont1_t' num2str(stepscale) '_' ascent '' s_w '' caseID '_' testnumber ];
                testcase0 = [ '_' num2str(N) '_LuCont0_t' num2str(stepscale) '_' ascent '' s_w '' caseID '_' testnumber ];

        end
        if sstest ~= 0 
            testcase_max = [ '_' num2str(N) '_L1Max_t' num2str(stepscale) '_tp(' num2str(tptest) ')_' ascent '' caseID ];
            testcase1 = [ '_' num2str(N) '_L1Cont1_t' num2str(stepscale) '_tp(' num2str(tptest) ')_' ascent '' caseID ]; 
            testcase0 = [ '_' num2str(N) '_L1Cont0_t' num2str(stepscale) '_tp(' num2str(tptest) ')_' ascent '' caseID ];
        end
        if ss_shift ~= 0
            testcase_max = [ '_' num2str(N) '_L1Max_t' num2str(stepscale) '_tp(' num2str(tptest) ')_shift' num2str(ss_shift) '_' ascent '' caseID ];
            testcase1 = [ '_' num2str(N) '_L1Cont1_t' num2str(stepscale) '_tp(' num2str(tptest) ')_shift' num2str(ss_shift) '_' ascent '' caseID ]; 
            testcase0 = [ '_' num2str(N) '_L1Cont0_t' num2str(stepscale) '_tp(' num2str(tptest) ')_shift' num2str(ss_shift) '_' ascent '' caseID ];
        end
        if shiftT ~= 0
            testcase_max = [ '_' num2str(N) '_LuMax_t' num2str(stepscale) '_shift' pm '' num2str(shiftT) '_' ascent '' caseID ];
            testcase1 = [ '_' num2str(N) '_LuCont1_t' num2str(stepscale) '_shift' pm '' num2str(shiftT) '_' ascent '' caseID ];
            testcase0 = [ '_' num2str(N) '_LuCont0_t' num2str(stepscale) '_shift' pm '' num2str(shiftT) '_' ascent '' caseID ];
        end
    else
        testcase_max = [ '_' num2str(N) '_slingshot' num2str(s) 'Max_t' num2str(stepscale) '_' ascent '' caseID ];
        testcase1 = [ '_' num2str(N) '_slingshot' num2str(s) 'Cont1_t' num2str(stepscale) '_' ascent '' caseID ];
        testcase0 = [ '_' num2str(N) '_slingshot' num2str(s) 'Cont0_t' num2str(stepscale) '_' ascent '' caseID ];
        if shiftT ~= 0
            testcase_max = [ '_' num2str(N) '_slingshot' num2str(s) 'Max_t' num2str(stepscale) ...
                '_shift' pm '' num2str(shiftT) '_' ascent '' caseID ];
            testcase1 = [ '_' num2str(N) '_slingshot' num2str(s) 'Cont1_t' num2str(stepscale) ...
                '_shift' pm '' num2str(shiftT) '_' ascent '' caseID ];
            testcase0 = [ '_' num2str(N) '_slingshot' num2str(s) 'Cont0_t' num2str(stepscale) ...
                '_shift' pm '' num2str(shiftT) '_' ascent '' caseID ];
        end
    end

return