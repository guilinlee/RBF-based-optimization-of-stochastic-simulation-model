function [hoptim, miu, inverseBRes]=RBFFit2(Xin,Ain,Rdatain)
%checked2
global X A D Res 
X=Xin; nX=size(X,1);
A=Ain; nA=size(A,1);
D=[repmat(X,nA,1) kron(A,ones(nX,1))];
Rdata=Rdatain;
options=optimset('Display','off','TolX',10^-12,'TolFun',10^-12);
startingpoint=[log(std(D,1,1)) log(var(Rdata)/20)];
miu=mean(Rdata);
Res=Rdata-miu;
[loghopt,fval,exitflag]=fminsearch(@crossvalidate,startingpoint,options);
hoptim=exp(loghopt);
nugget=hoptim(end);
[B1,B2]=RBF(hoptim(1:(end-1)));
%B1=E1L1E1^T; B2=E2L2E2^T
[E1, L1]=svd(B1);[E2, L2]=svd(B2);
%(B+nuggetI)^{-1}=(E1 Kronecker product kron E2)(nugget*I+L1 kron L2)(E1^T kron E2^T)
E1T=transpose(E1);E2T=transpose(E2);
L2L1I=diag(1./(diag(kron(L2,L1))+nugget));
inverseBRes=kron(E2,E1)*(L2L1I*((kron(E2T,E1T))*Res));

function CVE=crossvalidate(logh)
global D Res
h=exp(logh);
nugget=h(end);
if(sum(h(1:(end-1))>(12*std(D,1,1)))>=1)
    CVE=Inf;
    return
elseif(sum(h(1:(end-1))<(0.1*std(D,1,1)))>=1)
    CVE=Inf;
    return
end
%%to calculate (B+nuggetI)^{-1}, first calculate B1 B2 and 
[B1,B2]=RBF(h(1:(end-1)));
%B1=E1L1E1^T; B2=E2L2E2^T
[E1,L1]=svd(B1);[E2,L2]=svd(B2);
E1T=transpose(E1);E2T=transpose(E2);
L2L1I=diag(1./(diag(kron(L2,L1))+nugget));
KE2E1=kron(E2,E1);
DiaginverseB=sum((KE2E1').*(L2L1I*(kron(E2T,E1T))),1);
Error=diag(1./DiaginverseB)*(KE2E1*(L2L1I*((kron(E2T,E1T))*Res)));
CVE=mean((Error).^2);

function [B1,B2]=RBF(hin)
%%%%The RBF correlation matrix
global X A
n1=size(X,1);n2=size(A,1);
d1=size(X,2);d2=size(A,2);
h1=hin(1:d1);h2=hin((d1+1):(d1+d2));
B1=zeros(n1,n1);B2=zeros(n2,n2);
%%B1
for i=1:n1
    for j=i:n1
        z=(abs(X(i,:)-X(j,:)))./h1;
        B1(i,j)=prod(exp(-z).*(1+z),2);     
    end
end
for i=2:n1
    for j=1:(i-1)
        B1(i,j)=B1(j,i);
    end
end
%%B2
for i=1:n2
    for j=i:n2
        z=(abs(A(i,:)-A(j,:)))./h2;
        B2(i,j)=prod(exp(-z).*(1+z),2);     
    end
end
for i=2:n2
    for j=1:(i-1)
        B2(i,j)=B2(j,i);
    end
end