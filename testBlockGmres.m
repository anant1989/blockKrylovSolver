clear all;
A = dlmread('SystemMatrix_1.00e+006_real.txt') + 1i*dlmread('SystemMatrix_1.00e+006_imag.txt');
RHS1 = dlmread('System_1.00e+006_RHS_block1_real.txt') + 1i*dlmread('System_1.00e+006_RHS_block1_imag.txt') ;
RHS0 = dlmread('System_1.00e+006_RHS_block0_real.txt') + 1i*dlmread('System_1.00e+006_RHS_block0_imag.txt') ;

LHS0 = zeros(size(RHS0));
LHS1 = zeros(size(RHS0));
LHS0_ = dlmread('System_1.00e+006_LHS_block0_real.txt') + 1i*dlmread('System_1.00e+006_LHS_block0_imag.txt') ;
LHS1_ = dlmread('System_1.00e+006_LHS_block1_real.txt') + 1i*dlmread('System_1.00e+006_LHS_block1_imag.txt') ;

err = 0;
iter = 0;

%[LHS0,flag,relres,iter,resvec] = gmres(A, RHS0, length(RHS0), 1e-4);
%[LHS1,flag,relres,iter,resvec] = gmres(A, RHS1, length(RHS0), 1e-4);

%[LHS0, err, iter] = mygmres(A, LHS0, RHS0, length(RHS0), 1e-4);
%[LHS1, err, iter] = mygmres(A, LHS1, RHS1, length(RHS0), 1e-4);

%call block GMRES
LHS = [LHS0 LHS1];
RHS = [RHS0 RHS1];
%[LHS0, err, iter] = myblockgmres1(A, LHS(:,1), RHS(:,1), length(RHS0), 1e-4);
%[LHS1, err, iter] = myblockgmres1(A, LHS(:,2), RHS(:,2), length(RHS0), 1e-4);
[LHS, err, iter] = myblockgmres1(A, LHS, RHS, length(RHS0), 1e-4);
LHS0 = LHS(:,1);
LHS1 = LHS(:,2);

close all
plot(abs(LHS0_))
hold on
plot(abs(LHS0))
legend('expected LHS0','LHS0')
fprintf('relative error = %f\n', norm(LHS0-LHS0_)/norm(LHS0_))
fprintf('relative error = %f\n', norm(LHS1-LHS1_)/norm(LHS1_))
