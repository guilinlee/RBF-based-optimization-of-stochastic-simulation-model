function  [Qet] = release(u2,time,r0,u20,zL,nz,nr,dz,dr,r) 
%checked2
Qu1t=0;
%  for i=1:nz
%  for j=1:nr
%    Qu1t=Qu1t+u2(time,i,j)*r(j)*dr*dz;
%  end
%  end
  for i=1:(nz-1)
  for j=1:(nr-1)
    Qu1t=Qu1t+((u2(time,i,j)*r(j)+u2(time,i+1,j)*r(j)+u2(time,i,j+1)*r(j+1)+u2(time,i+1,j+1)*r(j+1))/4)*dr*dz;
  end
  end
% Qet=pi*r0^2*zL*u20-2*pi*Qu1t;
 Qet=pi*r0^2*zL*u20-2*2*pi*Qu1t;
