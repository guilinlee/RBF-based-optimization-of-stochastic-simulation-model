function [meanratio]=robustdesire(starvalue,optimdesignv, XRange, XMin,ARange, AMin, Design,hopt,miu,inverseBRes)
%checked2
x=starvalue;
sumratio=0;
L=size(optimdesignv,1);
for i=1:L
    A=optimdesignv(i,3:end);
    desireA=optimdesignv(i,2);
    desirestarA=-nRBFPredict((x-XMin)./XRange,(A-AMin)./ARange,Design,hopt,miu,inverseBRes);
    ratioA=desirestarA/desireA;
    sumratio=sumratio+ratioA;
end
meanratio=sumratio/L;
