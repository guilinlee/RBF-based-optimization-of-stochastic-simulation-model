function arf=computermodel(designr0)
% the output of the simulator is the F1 in the loss function 
%%%%%%%%%%%%%%%%%%%checked2
% Dynamic analysis of drug delivery
%
% The following equations model a polymer matrix in cylindrical
% coordinates 
%
%   u1 material balance
%
%   u1t = Du1*(u1rr + (1/r)u1r)                                 (1)
%
%   u2 material balance
%
%   u2t = Du2*(u2rr + (1/r)*u2r)                                (2)
%
% The variables and parameters for this model are (in cgs units)
%
%   Water concentration            u1                     (eq. (1))
%
%   Drug concentration             u2                     (eq. (2))
%
%   Time                            t
%
%   Radial position                 r
%
%   Axial position                  z
%
%   Polymer matrix radius          r0              1
%
%   Polymer matrix length          zL              2
%
%   Initial u1                    u10            0.5
%
%   Initial u2                    u20              1
%
%   External u1                   u1e              1
%
%   External u2                   u2e              0
%
%   u1 diffusivity                Du1        1.0e-06
%
%   u2 diffusivity                Du2        1.0e-06
%
%   u1 mass transfer coefficient  ku1            0.1
%
%   u2 mass transfer coefficient  ku2            0.1
%
% The solution to this system, u1(r,z,t) (from eq. (1)) 
% and u2(r,z,t) (from eq. (2)) is computed by ode15s for 
% the integration in t. 
% Global area
%
% Global area
  global     nr     nz     dr     dz    drs    dzs...
              r     r0      z     zL    Du1    Du2...
             u1     u2    u1e    u2e    ku1    ku2...
          ncall
%%%%%%%%%
%desired release profile

% Model parameters
%random parameter the concentration of drug in the environment is random
%uniform distribution from [0,0.5]
  %u2e=0;
u2e=unifrnd(0,0.5);
%design parameters
 r0=designr0; zL=2;
 %%r0=1;zL=2
  u10=0.5; u20=1;u1e=1;
  Du1=1.0e-06; Du2=1.0e-06;
  ku1=1.0e-01; ku2=unifrnd(0.05,0.15);
%
% Grid in axial direction
  nz=11;
  dz=zL/(2*(nz-1));
  for i=1:nz
    z(i)=zL/2+(i-1)*dz;
  end
  dzs=dz^2;
%
% Grid in radial direction
  nr=11;
  dr=r0/(nr-1);
  for j=1:nr
    r(j)=(j-1)*dr;
  end
  drs=dr^2;  
%
% Independent variable for ODE integration
  tf=2*3600*24;
  tout=[0.0:2*900*24:tf]'; 
  nout=5;
  ncall=0;
%
% Initial condition
  for i=1:nz
  for j=1:nr
    u1(i,j)=u10;
    u2(i,j)=u20;
    u0((i-1)*nr+j)=u1(i,j);
    u0((i-1)*nr+j+nz*nr)=u2(i,j);
  end
  end
%
% ODE integration 
  reltol=1.0e-04; abstol=1.0e-04;
  options=odeset('RelTol',reltol,'AbsTol',abstol);
  [t,u]=ode15s(@pde_1,tout,u0,options);
%
% 1D to 2D matrices
  for it=1:nout
  for i=1:nz
  for j=1:nr  
    u1(it,i,j)=u(it,(i-1)*nr+j);
    u2(it,i,j)=u(it,(i-1)*nr+j+nz*nr);
  end
  end
  end

  %%%%%%%%%%%%the total amount of the drug that has left the PM at time t
  %%%%%%%%%%%%Q_e(t) is 
  %%desired drug release profile
%%
  arf=zeros(5,1);
  for i=1:5
  time=i;
  arf(i)=release(u2,time,r0,u20,zL,nz,nr,dz,dr,r);
  end
  %%
end