%function Example
%% deterministic design input: r0
%checked2
XMax=1.5;
XMin=0.5;
XRange=XMax-XMin;
XDesign=[0:(1/4):1]';
n1=size(XDesign,1);
UncodedXDesign=XDesign*XRange+XMin;
AMax=3;
AMin=0.5;
ARange=AMax-AMin;
ADesign=[0:(1/20):1]';
n2=size(ADesign,1);
UncodedADesign=ADesign*ARange+AMin;
Design=[repmat(XDesign,n2,1) kron(ADesign,ones(n1,1))];
UncodedDesign=[repmat(UncodedXDesign,n2,1) kron(UncodedADesign,ones(n1,1))];
n=size(Design,1);
%% number of design point
%%%%%number of observations at each design point
Rep=5;
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
UncodedDesign2=[kron(UncodedXDesign,ones(n2,1)) repmat(UncodedADesign,n1,1)];
figure(100),plot3(UncodedDesign2(:,1),UncodedDesign2(:,2),Ldata,'o','color','red'),hold on
xlabel('$x$','interpreter','latex'),ylabel('$\theta_1$','interpreter','latex'),zlabel('Desirability score')
for i=1:n
    figure(100),plot3(repmat(UncodedDesign2(i,1),1,2),repmat(UncodedDesign2(i,2),1,2),[min(Ldata(i,:)) max(Ldata(i,:))],'LineWidth',1,'color','blue')
end
%%RBF fit to the Rdata from simulator  [hopt, miu, RIRes]=RBFFit2(Xin,Ain,Rin)
[hopt, miu, inverseBRes]=RBFFit2(XDesign,ADesign,Rdata);
RhatDesign=RBFPredict2(Design,Design,hopt,miu,inverseBRes);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prediction of testset
Xtest=[0:(1/20):1]';Atest=[0:(1/20):1]';
UncodedXtest=Xtest*XRange+XMin;
UncodedAtest=Atest*ARange+AMin;
ntestX=size(Xtest,1);ntestA=size(Atest,1);
TestSet=[repmat(Xtest,ntestA,1) kron(Atest,ones(ntestX,1))];
UncodedTestSet=[repmat(UncodedXtest,ntestA,1) kron(UncodedAtest,ones(ntestX,1))];
ntest=size(TestSet,1);
RhatTest=RBFPredict2(TestSet,Design,hopt,miu,inverseBRes);
figure(1),mesh(UncodedXtest,UncodedAtest,reshape(RhatTest,ntestX,ntestA)')
xlabel('$x$','interpreter','latex'),ylabel('$\theta_1$','interpreter','latex'),zlabel('Desirability score')

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
setA=[(AMin+0.01/2):0.01:(AMax-0.01/2)]';
LA=size(setA,1);
optimdesignfit=zeros(LA,3);
for i=1:LA
f = @(x)nRBFPredict((x-XMin)/XRange,(setA(i)-AMin)/ARange,Design,hopt,miu,inverseBRes);
[optimdesign,optimdesire]=patternsearch(f,(XMin+XMax)/2,[],[],[],[],XMin,XMax);
optimdesignfit(i,:)=[optimdesign,-optimdesire,setA(i)];
end
optimdesignfit
figure(2),plot(setA,optimdesignfit(:,1),'LineWidth',2)
xlabel('$\theta_1$','interpreter','latex'),ylabel('$x^*(\theta_1)$','interpreter','latex')
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
figure(3),plot(optimdesignfit(:,3), meanratiolist)
xlabel('$\theta_1$','interpreter','latex'),ylabel('$Q(\theta_1)$','interpreter','latex')
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
figure(4), hold on
h=plot3d_errorbars(UncodedTestSet(:,1), UncodedTestSet(:,2), RhatTest, UCL, LCL, RtrueTest)
xlabel('$x$','interpreter','latex'),ylabel('$\theta_1$','interpreter','latex')

indices=[];
for i=1:ntestA
    if(mod(i,2)==1)
        indices=[indices ((i-1)*ntestX+1):(i*ntestX)];
    else
        indices=[indices (i*ntestX):-1:((i-1)*ntestX+1)];
    end
end
figure(5), plot(1:ntest,RhatTest(indices),'LineWidth',1,'color','red'), hold on, 
plot(1:ntest,RtrueTest(indices),'.','MarkerSize',10,'color','black')
plot(1:ntest,LCL(indices),':','LineWidth',1,'color','blue')
plot(1:ntest,UCL(indices),':','LineWidth',1,'color','blue')
legend('Predicted desirability score','Accurate estimate of desirability score','LPIL','UPIL','location','best')
xlabel('Grid point label'),ylabel('Desirability score')
