function [invR, logdetR]=invandlogdet(R)
%checked 2
%
[CR, p]=chol(R);
if(p~=0)
    [Eigvec, Eigval]=svd(R);
    Eigval=diag(Eigval);
    invR=Eigvec*diag(1./Eigval)*Eigvec';
    logdetR=sum(log(Eigval));
else
    ICR=inv(CR);
    invR=ICR*(ICR');
    logdetR=2*sum(log(diag(CR)));
end