  function ut=pde_3(t,u)
%
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
% Variable diffusivities
  ts=24*60*60;
  Du1t=Du1*(1+t/ts);
  Du2t=Du2*(1+t/ts);
%
% Step through the grid points in r
  for i=1:nz
    u1_1d=u1(i,:);  
    u2_1d=u2(i,:);      
%
%   u1r, u2r
    u1r_1d=dss004(0,r0,nr,u1_1d);
    u1r(i,:)=u1r_1d;
    u1r(i,1)= 0.0;
    u1r(i,nr)=(ku1/Du1t)*(u1e-u1(i,nr));    
    u2r_1d=dss004(0,r0,nr,u2_1d);
    u2r(i,:)=u2r_1d;
    u2r(i,1)= 0;
    u2r(i,nr)=(ku2/Du2t)*(u2e-u2(i,nr));
%
%   u1rr, u2rr
    u1r_1d( 1)=0;
    u1r_1d(nr)=(ku1/Du1t)*(u1e-u1_1d(nr));
    nl=2; nu=2;
    u1rr_1d=dss044(0,r0,nr,u1_1d,u1r_1d,nl,nu);
    u1rr(i,:)=u1rr_1d;
    u2r_1d( 1)=0;
    u2r_1d(nr)=(ku2/Du2t)*(u2e-u2_1d(nr));
    nl=2; nu=2;
    u2rr_1d=dss044(0,r0,nr,u2_1d,u2r_1d,nl,nu);
    u2rr(i,:)=u2rr_1d; 
%
%   (1/r)*u1r, (1/r)*u2r    
    for j=1:nr 
      if(j~=1)
        u1r(i,j)=(1.0/r(j))*u1r(i,j);
        u2r(i,j)=(1.0/r(j))*u2r(i,j);
      elseif(j==1) 
        u1rr(i,j)=2.0*u1rr(i,j);
        u2rr(i,j)=2.0*u2rr(i,j);
      end 
%
%   Next j
    end
%
% Next i
  end
%
% Step through the grid points in z
  for j=1:nr
    u1_1d=u1(:,j);  
    u2_1d=u2(:,j);      
%
%   u1zz, u2zz
    u1z_1d( 1)=0.0;
    u1z_1d(nz)=(ku1/Du1t)*(u1e-u1_1d(nz));
    nl=2; nu=2;
    u1zz_1d=dss044(zL/2,zL,nz,u1_1d,u1z_1d,nl,nu);
    u1zz(:,j)=u1zz_1d;
    u2z_1d( 1)=0.0;
    u2z_1d(nz)=(ku2/Du2t)*(u2e-u2_1d(nz));
    nl=2; nu=2;
    u2zz_1d=dss044(zL/2,zL,nz,u2_1d,u2z_1d,nl,nu);
    u2zz(:,j)=u2zz_1d; 
%
% Next j
  end
%
% PDEs
  for i=1:nz
  for j=1:nr
      u1t(i,j)=Du1t*(u1rr(i,j)+u1r(i,j)+u1zz(i,j));
      u2t(i,j)=Du2t*(u2rr(i,j)+u2r(i,j)+u2zz(i,j));
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