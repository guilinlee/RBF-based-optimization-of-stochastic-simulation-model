function d=dMin(fx,lower,upper,A)
%checked2
if(fx>upper)
    d=0;
elseif((fx>=lower) && (fx<=upper))
    d=((fx-upper)/(lower-upper))^A;
else
    d=1;
end