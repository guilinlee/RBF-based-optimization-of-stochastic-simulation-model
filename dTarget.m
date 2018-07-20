function d=dTarget(fx,lower,target,upper)
%checked2
if((fx>=lower) && (fx<=target))
    d=((fx-lower)/(target-lower))^2;
elseif((fx>=target) && (fx<=upper))
    d=((fx-upper)/(target-upper));
else
    d=0;
end