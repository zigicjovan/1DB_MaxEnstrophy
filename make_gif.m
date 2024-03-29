E0 = 1000;
TW = 0.028563;

filename = [ 'optICphys_4096_slingshot1_shiftp00_E0(' num2str(E0) ')_Timept_14_lambda(0.5)' ];
casename = 'slingshot1';
solution_in_time_figure(E0,filename,TW,casename)

filename = [ 'optICphys_4096_slingshot2_shiftp00_E0(' num2str(E0) ')_Timept_14_lambda(0.5)' ];
casename = 'slingshot2';
solution_in_time_figure(E0,filename,TW,casename)

E0 = 100000;
TW = 0.0026522;

filename = [ 'optICphys_4096_Lu_t4_E0(' num2str(E0) ')_Timept_13_lambda(0.5)' ];
casename = 'original';
solution_in_time_figure(E0,filename,TW,casename)

filename = [ 'optICphys_4096_slingshot1_t10_E0(' num2str(E0) ')_Timept_13_lambda(0.5)' ];
casename = 'slingshot1';
solution_in_time_figure(E0,filename,TW,casename)

filename = [ 'optICphys_4096_slingshot2_t10_E0(' num2str(E0) ')_Timept_13_lambda(0.5)' ];
casename = 'slingshot2';
solution_in_time_figure(E0,filename,TW,casename)

filename = [ 'optICphys_4096_slingshot3_t2_E0(' num2str(E0) ')_Timept_13_lambda(0.5)' ];
casename = 'slingshot3';
solution_in_time_figure(E0,filename,TW,casename)

filename = [ 'optICphys_4096_slingshot4_t2_E0(' num2str(E0) ')_Timept_13_lambda(0.5)' ];
casename = 'slingshot4';
solution_in_time_figure(E0,filename,TW,casename)

filename = [ 'optICphys_4096_slingshot5_t10_E0(' num2str(E0) ')_Timept_13_lambda(0.5)' ];
casename = 'slingshot5a';
solution_in_time_figure(E0,filename,TW,casename)

filename = [ 'optICphys_4096_slingshot6_t4_E0(' num2str(E0) ')_Timept_13_lambda(0.5)' ];
casename = 'slingshot5b';
solution_in_time_figure(E0,filename,TW,casename)