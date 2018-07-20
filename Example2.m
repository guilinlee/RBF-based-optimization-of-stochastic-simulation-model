%function Example2
%% deterministic design input: r0
%checked2
XMax=1.5;
XMin=0.5;
XRange=XMax-XMin;
XDesign=[0:(1/4):1]';
n1=size(XDesign,1);
UncodedXDesign=XDesign*XRange+XMin;
AMax=[3 1.25];
AMin=[0.5 0.75];
ARange=AMax-AMin;
ADesign=[repmat([0:(1/10):1]',6,1) kron([0:(1/5):1]',ones(11,1))];
nX=size(XDesign,1); nA=size(ADesign,1);
UncodedADesign=ADesign.*repmat(ARange,nA,1)+repmat(AMin,nA,1);
Design=[repmat(XDesign,nA,1) kron(ADesign,ones(nX,1))];
%% number of design point
%%%%%number of observations at each design point
Rep=10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% generate data Y|X from the computer model
Data=cell(1,n1);
for i=1:n1
    Data{i}=zeros(Rep,5);
end
for i=1:n1
    for repi=1:Rep
         Data{i}(repi,:)=computermodel((UncodedXDesign(i,:)));
    end
end
[Rdata,Ldata]=computeloss(Data,UncodedXDesign,UncodedADesign);
%%RBF fit to the Rdata from simulator  [hopt, miu, RIRes]=RBFFit2(Xin,Ain,Rin)
[hopt, miu, inverseBRes]=RBFFit2(XDesign,ADesign,Rdata);
%RhatDesign=RBFPredict2(Design,Design,hopt,miu,inverseBRes);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prediction of testset
Xtest=[0:(1/20):1]';Atest=[repmat([0:(1/10):1]',11,1) kron([0:(1/10):1]',ones(11,1))];
ntestX=size(Xtest,1);ntestA=size(Atest,1);
UncodedXtest=Xtest*XRange+XMin;
UncodedAtest=Atest.*repmat(ARange,ntestA,1)+repmat(AMin,ntestA,1);
TestSet=[repmat(Xtest,ntestA,1) kron(Atest,ones(ntestX,1))];
ntest=size(TestSet,1);
RhatTest=RBFPredict2(TestSet,Design,hopt,miu,inverseBRes);

RepTest=100;
DataTest=cell(1,ntestX);
for i=1:ntestX
    DataTest{i}=zeros(RepTest,5);
end
for i=1:ntestX
    for repi=1:RepTest
         DataTest{i}(repi,:)=computermodel((UncodedXtest(i,:)));
    end
end
RtrueTest=computeloss(DataTest,UncodedXtest,UncodedAtest);

MSE=mean((RtrueTest-RhatTest).^2)
%%%%%
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%optimal design point for different A
da=length(AMin); nostepsA=20;
setA=(fullfact(nostepsA*ones(1,da))-0.5)/nostepsA; LA=nostepsA^da;
setA=setA.*repmat(ARange,LA,1)+repmat(AMin,LA,1);
optimdesignfit=zeros(LA,4);
for i=1:LA
f = @(x)nRBFPredict((x-XMin)/XRange,(setA(i,:)-AMin)./ARange,Design,hopt,miu,inverseBRes);
[optimdesign,optimdesire]=patternsearch(f,(XMin+XMax)/2,[],[],[],[],XMin,XMax);
optimdesignfit(i,:)=[optimdesign,-optimdesire,setA(i,:)];
end
optimdesignfit
figure(2),mesh(unique(setA(:,1)),unique(setA(:,2)),reshape(optimdesignfit(:,1),nostepsA,nostepsA)')
xlabel('$\theta_1$','interpreter','latex'),ylabel('$\theta_2$','interpreter','latex'),zlabel('$x^*(\theta)$','interpreter','latex')
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
%%calculate the robustness measure for each optimal solution
meanratiolist=zeros(LA,1);
for j=1:LA
 starvalue=optimdesignfit(j,1);
 [meanratio]=robustdesire(starvalue, optimdesignfit, XRange, XMin,ARange, AMin, Design,hopt,miu,inverseBRes);
 meanratiolist(j)=meanratio;
end
%%%
figure(3),mesh(unique(optimdesignfit(:,3)),unique(optimdesignfit(:,4)),reshape(meanratiolist,nostepsA,nostepsA)')
xlabel('$\theta_1$','interpreter','latex'),ylabel('$\theta_2$','interpreter','latex'),zlabel('$Q(\theta)$','interpreter','latex')
r=find(meanratiolist==max(meanratiolist));
optimdesignfit(r,:)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%bootstrap to get the expected loss R 
nboot=1000;
ExpectedRBoot=zeros(ntest,nboot);
BootData=cell(1,n1);
for i=1:n1
    BootData{i}=zeros(Rep,5);
end
for i=1:nboot
    for j=1:n1
        indices=randsample(1:Rep,Rep,'true');
        BootData{j}=Data{j}(indices,:);
    end
    Rdataboot=computeloss(BootData,UncodedXDesign,UncodedADesign);
    [bhopt, bmiu, binverseBRes]=RBFFit2(XDesign,ADesign,Rdataboot);
    bRhat=RBFPredict2(TestSet,Design,bhopt,bmiu,binverseBRes);
    ExpectedRBoot(:,i)=bRhat;
end

MeanExpectedRBoot=mean(ExpectedRBoot,2);
UCL=prctile(ExpectedRBoot,99,2);
LCL=prctile(ExpectedRBoot,1,2);

%%%%%%
coverage=mean((LCL<=RtrueTest).*(RtrueTest<=UCL))
indices=[];
for i=1:11
    for j=1:11
    if(mod(i,2)==1)
        if(mod(j,2)==1)
            indices=[indices ((i-1)*11*ntestX+(j-1)*ntestX+1):((i-1)*11*ntestX+j*ntestX)];
        else
            indices=[indices ((i-1)*11*ntestX+j*ntestX):-1:((i-1)*11*ntestX+(j-1)*ntestX+1)];
        end
    else
        if(mod(j,2)==1)
            indices=[indices (i*11*ntestX-(j-1)*ntestX):-1:(i*11*ntestX-j*ntestX+1)];
        else
            indices=[indices (i*11*ntestX-j*ntestX+1):(i*11*ntestX-(j-1)*ntestX)];
        end       
    end
    end
end
figure(5), plot(1:ntest,RhatTest(indices),'LineWidth',1,'color','red'), hold on, 
plot(1:ntest,RtrueTest(indices),'.','MarkerSize',10,'color','black')
plot(1:ntest,LCL(indices),':','LineWidth',1,'color','blue')
plot(1:ntest,UCL(indices),':','LineWidth',1,'color','blue')
legend('Predicted desirability score','Accurate estimate of desirability score','LCL','UCL','location','best')
xlabel('Grid point label'),ylabel('Desirability score')
