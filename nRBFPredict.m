function [nyhat]=nRBFPredict(x,A,D,hopt,miu,inverseBRes)
%checked
TestSet=[x(:) repmat(A,length(x),1)];
[yhat]=RBFPredict2(TestSet,D,hopt,miu,inverseBRes);
nyhat=-yhat;
