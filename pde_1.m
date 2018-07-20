  function ut=pde_1(t,u)
%checked2
% Global area
  global     nr     nz     dr     dz    drs    dzs...
              r     r0      z     zL    Du1    Du2...
             u1     u2    u1e    u2e     ku1   ku2...
          ncall
%
% 1D to 2D matrices
  for i=1:nz
  for j=1:nr
    ij=(i-1)*nr+j;
    u1(i,j)=u(ij);
    u2(i,j)=u(ij+nr*nz);
  end
  end
%
% Step through the grid points in r and z
  for i=1:nz
  for j=1:nr  
%
%   (1/r)*u1r, (1/r)*u2r
    if(j==1)
      u1r(i,j)=2*(u1(i,j+1)-u1(i,j))/drs;
      u2r(i,j)=2*(u2(i,j+1)-u2(i,j))/drs;
    elseif(j==nr)
      u1r(i,j)=(1/r(j))*(ku1/Du1)*(u1e-u1(i,j));
      u2r(i,j)=(1/r(j))*(ku2/Du2)*(u2e-u2(i,j));
    else
      u1r(i,j)=(1/r(j))*(u1(i,j+1)-u1(i,j-1))/(2*dr);
      u2r(i,j)=(1/r(j))*(u2(i,j+1)-u2(i,j-1))/(2*dr);
    end      
%
%   u1rr, u2rr
    if(j==1)
      u1rr(i,j)=2*(u1(i,j+1)-u1(i,j))/drs;
      u2rr(i,j)=2*(u2(i,j+1)-u2(i,j))/drs;
    elseif(j==nr)
      u1f=u1(i,j-1)+2*dr*ku1/Du1*(u1e-u1(i,j));
      u1rr(i,j)=(u1f-2*u1(i,j)+u1(i,j-1))/drs;
      u2f=u2(i,j-1)+2*dr*ku2/Du2*(u2e-u2(i,j));
      u2rr(i,j)=(u2f-2*u2(i,j)+u2(i,j-1))/drs;
    else
      u1rr(i,j)=(u1(i,j+1)-2*u1(i,j)+u1(i,j-1))/drs;
      u2rr(i,j)=(u2(i,j+1)-2*u2(i,j)+u2(i,j-1))/drs;
    end
%
%   u1zz, u2zz
    if(i==1)
      u1zz(i,j)=2*(u1(i+1,j)-u1(i,j))/dzs;
      u2zz(i,j)=2*(u2(i+1,j)-u2(i,j))/dzs;        
    elseif(i==nz)
      u1f=u1(i-1,j)+2*dz*ku1/Du1*(u1e-u1(i,j));
      u1zz(i,j)=(u1f-2*u1(i,j)+u1(i-1,j))/dzs;
      u2f=u2(i-1,j)+2*dz*ku2/Du2*(u2e-u2(i,j));
      u2zz(i,j)=(u2f-2*u2(i,j)+u2(i-1,j))/dzs;
    else
      u1zz(i,j)=(u1(i+1,j)-2*u1(i,j)+u1(i-1,j))/dzs;
      u2zz(i,j)=(u2(i+1,j)-2*u2(i,j)+u2(i-1,j))/dzs;
    end    
%
%   PDEs
      u1t(i,j)=Du1*(u1rr(i,j)+u1r(i,j)+u1zz(i,j));
      u2t(i,j)=Du2*(u2rr(i,j)+u2r(i,j)+u2zz(i,j));
  end
  end
%
% 2D to 1D matrices
  for i=1:nz
  for j=1:nr
    ij=(i-1)*nr+j;
    ut(ij)=u1t(i,j);
    ut(ij+nr*nz)=u2t(i,j);
  end
  end 
%
% Transpose and count
  ut=ut';
  ncall=ncall+1;