E0 = 1000;
% TW = 0.0026522;
TW = 0.026522;

filename0 = [ 'optICphys_2048_Lu_t4_E0(' num2str(E0) ')_Timept_13_lambda(0.5)' ];
filename1 = [ 'optICphys_4096_slingshot1_t2_E0(' num2str(E0) ')_Timept_13_lambda(0.5)' ];
filename2 = [ 'optICphys_4096_slingshot2_t2_E0(' num2str(E0) ')_Timept_13_lambda(0.5)' ];
filename3 = [ 'optICphys_4096_slingshot3_t2_E0(' num2str(E0) ')_Timept_13_lambda(0.5)' ];
filename4 = [ 'optICphys_4096_slingshot4_t2_E0(' num2str(E0) ')_Timept_13_lambda(0.5)' ];
filename5 = [ 'optICphys_4096_slingshot5_t2_E0(' num2str(E0) ')_Timept_13_lambda(0.5)' ];
filename6 = [ 'optICphys_4096_slingshot6_t2_E0(' num2str(E0) ')_Timept_13_lambda(0.5)' ];

solution_derivative(E0,filename0,TW);
solution_derivative(E0,filename1,TW);
solution_derivative(E0,filename2,TW);
solution_derivative(E0,filename3,TW);
solution_derivative(E0,filename4,TW);
solution_derivative(E0,filename5,TW);
solution_derivative(E0,filename6,TW);