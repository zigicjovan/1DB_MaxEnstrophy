function [testcase_max,testcase1,testcase0] = testcase_name(N,stepscale,ascent,testnumber)   %
                                                                                                                               %                           %
    testcase_max = [ '_' num2str(N) '_LuMax_t' num2str(stepscale) '_' ascent '' num2str(testnumber) ];       %
    testcase1 = [ '_' num2str(N) '_LuCont1_t' num2str(stepscale) '_' ascent '' num2str(testnumber) ];        %
    testcase0 = [ '_' num2str(N) '_LuCont0_t' num2str(stepscale) '_' ascent '' num2str(testnumber) ];        %                                                                                                                 %

return