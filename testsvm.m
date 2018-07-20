
nsim=10;
SVRMAE=zeros(nsim,1);
SVRMSE=zeros(nsim,1);
time1=zeros(nsim,1);

for i=1:nsim
tic
%Mdl = fitrsvm(Design,Rdata,'standardize',false,'OptimizeHyperparameters',{'BoxConstraint','KernelScale','Epsilon','PolynomialOrder'},'KernelFunction','polynomial','HyperparameterOptimizationOptions',struct('Optimizer','gridsearch'))
%Mdl = fitrsvm(Design,Rdata,'standardize',false,'OptimizeHyperparameters',{'BoxConstraint','KernelScale','Epsilon'},'KernelFunction','polynomial','HyperparameterOptimizationOptions',struct('Optimizer','gridsearch'))
%Mdl = fitrsvm(Design,Rdata,'standardize',false,'OptimizeHyperparameters',{'BoxConstraint','KernelScale','Epsilon','PolynomialOrder'},'KernelFunction','polynomial')
Mdl = fitrsvm(Design,Rdata,'standardize',false,'OptimizeHyperparameters',{'BoxConstraint','KernelScale','Epsilon'},'KernelFunction','polynomial')
%Mdl.HyperparameterOptimizationResults.XTrace;

SVRypredtest = predict(Mdl,TestSet);

time1(i)=toc

SVRPredErr=RtrueTest-SVRypredtest;

SVRMAE(i)=mean(abs(SVRPredErr));
SVRMSE(i)=mean(SVRPredErr.^2);

end

return

tic
[hopt, miu, inverseBRes]=RBFFit2(XDesign,ADesign,Rdata);
RhatTest=RBFPredict2(TestSet,Design,hopt,miu,inverseBRes);
MSE=mean((RtrueTest-RhatTest).^2)
MAE=mean(abs(RtrueTest-RhatTest))
time0=toc

return
%Mdl2 = fitrsvm(Design,Rdata,'standardize',false,'OptimizeHyperparameters',{'BoxConstraint','KernelScale','Epsilon'},'KernelFunction','linear',struct('Optimizer','gridsearch'))
Mdl2 = fitrsvm(Design,Rdata,'standardize',false,'OptimizeHyperparameters',{'BoxConstraint','KernelScale','Epsilon'},'KernelFunction','linear','HyperparameterOptimizationOptions',struct('MaxObjectiveEvaluations',100))

%Mdl2.HyperparameterOptimizationResults.XTrace;

SVRypredtest2 = predict(Mdl2,TestSet);

time2=toc

SVRPredErr2=RtrueTest-SVRypredtest2;

SVRMAE2=mean(abs(SVRPredErr2))
SVRMSE2=mean(SVRPredErr2.^2)

