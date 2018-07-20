function [yhat]=RBFPredict2(TestSet,D,hopt,miu,inverseBRes)
%checked2
r=RBF2(TestSet,D,hopt);
yhat=miu+r*inverseBRes;

function Matrix=RBF2(x1,x2,hin)
h=hin(1:(end-1));
nugget=hin(end);
m1=size(x1,1);
m2=size(x2,1);
Matrix=zeros(m1,m2);
for i=1:m1
    for j=1:m2
        z=(abs(x1(i,:)-x2(j,:)))./h;
        Matrix(i,j)=prod(exp(-z).*(1+z),2); 
    end
end