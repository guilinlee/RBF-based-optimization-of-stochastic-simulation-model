function [Rdata,Ldata]=computeloss(Data,UncodedX,UncodedA)
%checked2
X=UncodedX;A=UncodedA;
n1=size(X,1);n2=size(A,1);
Ldata=[];
for i=1:n1
    for j=1:n2
        Ldata=[Ldata; lossfun(Data{i},X(i,:),A(j,:))'];
        Rdata(i,j)=mean(Ldata(end,:));
    end
end
Rdata=Rdata(:);

function loss=lossfun(Data,X,A)
if(length(A)==1)
    d24para=[1,3.5,12]; 
else
    d24para=[A(2),3.5,12];     
end
d48para=[1.2,7,13];
r0=X;
u20=1;
zL=2;
drug=pi*r0^2*zL*u20;
F2=dMin(drug,0.1,15,A(1));

Rep=size(Data,1);
loss=zeros(Rep,1);
for i=1:Rep
d24=dTarget(Data(i,3),d24para(1),d24para(2),d24para(3));
d48=dTarget(Data(i,5),d48para(1),d48para(2),d48para(3));
F1=d24*d48;
loss(i)=(F1.*F2).^(1/3);
end