  function ut=pde_4(t,u)
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
% Step through the grid points in r and z
  for i=1:nz
  for j=1:nr  
%
%   First RHS term in r
%
%   Constant Du1, Du2
%   Du1=Du1; Du2=Du2; dDu1=0; dDu2=0; 
%
%   Variable Du1, Du2
    u1crit=1; beta_u1=1; Du1crit=1.0e-06;
              beta_u2=1; Du2crit=1.0e-06;
    Du1=Du1crit*exp(-beta_u1*(1-u1(i,j)/u1crit)); 
    Du2=Du2crit*exp(-beta_u2*(1-u1(i,j)/u1crit));    
    dDu1=Du1*(beta_u1/u1crit); dDu2=0;     
    if(j==1)
      term1_1(i,j)=Du1*2*(u1(i,j+1)-u1(i,j))/drs;
      term1_2(i,j)=Du2*2*(u2(i,j+1)-u2(i,j))/drs;
    elseif(j==nr)
      u1f=u1(i,j-1)+2*dr*(ku1/Du1)*(u1e-u1(i,j));  
      term1_1(i,j)=Du1*(u1f-2*u1(i,j)+u1(i,j-1))/drs;
      u2f=u2(i,j-1)+2*dr*(ku2/Du2)*(u2e-u2(i,j)); 
      term1_2(i,j)=Du2*(u2f-2*u2(i,j)+u2(i,j-1))/drs;      
    else
      term1_1(i,j)=Du1*(u1(i,j+1)-2*u1(i,j)+u1(i,j-1))/drs;
      term1_2(i,j)=Du2*(u2(i,j+1)-2*u2(i,j)+u2(i,j-1))/drs;
    end 
%
%   Second RHS term in r
    if(j==1)
      term2_1(i,j)=0;
      term2_2(i,j)=0;
    elseif(j==nr)
      u1f=u1(i,j-1)+2*dr*(ku1/Du1)*(u1e-u1(i,j));  
      term2_1(i,j)=dDu1*((u1f-u1(i,j-1))/(2*dr))^2;
      u2f=u2(i,j-1)+2*dr*(ku2/Du2)*(u2e-u2(i,j)); 
      term2_2(i,j)=dDu2*((u2f-u2(i,j-1))/(2*dr))^2;      
    else
      term2_1(i,j)=dDu1*((u1(i,j+1)-u1(i,j-1))/(2*dr))^2;
      term2_2(i,j)=dDu2*((u2(i,j+1)-u2(i,j-1))/(2*dr))^2;
    end  
%
%   Third RHS term in r
    if(j==1)
      term3_1(i,j)=term1_1(i,j);
      term3_2(i,j)=term1_2(i,j);
    elseif(j==nr)
      u1f=u1(i,j-1)+2*dr*(ku1/Du1)*(u1e-u1(i,j));  
      term3_1(i,j)=(Du1/r(j))*(u1f-u1(i,j-1))/(2*dr);
      u2f=u2(i,j-1)+2*dr*(ku2/Du2)*(u2e-u2(i,j)); 
      term3_2(i,j)=(Du2/r(j))*(u2f-u2(i,j-1))/(2*dr);      
    else
      term3_1(i,j)=(Du1/r(j))*(u1(i,j+1)-u1(i,j-1))/(2*dr);
      term3_2(i,j)=(Du2/r(j))*(u2(i,j+1)-u2(i,j-1))/(2*dr);
    end          
%
%   First RHS term in z
    if(i==1)
      term4_1(i,j)=Du1*2*(u1(i+1,j)-u1(i,j))/dzs;
      term4_2(i,j)=Du2*2*(u2(i+1,j)-u2(i,j))/dzs;
    elseif(i==nz)
      u1f=u1(i-1,j)+2*dz*(ku1/Du1)*(u1e-u1(i,j));  
      term4_1(i,j)=Du1*(u1f-2*u1(i,j)+u1(i-1,j))/dzs;
      u2f=u2(i-1,j)+2*dz*(ku2/Du2)*(u2e-u2(i,j)); 
      term4_2(i,j)=Du2*(u2f-2*u2(i,j)+u2(i-1,j))/dzs;      
    else
      term4_1(i,j)=Du1*(u1(i+1,j)-2*u1(i,j)+u1(i-1,j))/dzs;
      term4_2(i,j)=Du2*(u2(i+1,j)-2*u2(i,j)+u2(i-1,j))/dzs;
    end 
%
%   Second RHS term in z
    if(i==1)
      term5_1(i,j)=0;
      term5_2(i,j)=0;
    elseif(i==nz)
      u1f=u1(i-1,j)+2*dz*(ku1/Du1)*(u1e-u1(i,j));  
      term5_1(i,j)=dDu1*((u1f-u1(i-1,j))/(2*dz))^2;
      u2f=u2(i-1,j)+2*dz*(ku2/Du2)*(u2e-u2(i,j)); 
      term5_2(i,j)=dDu2*((u2f-u2(i-1,j))/(2*dz))^2;      
    else
      term5_1(i,j)=dDu1*((u1(i+1,j)-u1(i-1,j))/(2*dz))^2;
      term5_2(i,j)=dDu2*((u2(i+1,j)-u2(i-1,j))/(2*dz))^2;
    end 
%
% Next j
  end
%
% Next i
  end
%
%   PDEs
  for i=1:nz
  for j=1:nr
      u1t(i,j)=term1_1(i,j)+term2_1(i,j)+term3_1(i,j)...
              +term4_1(i,j)+term5_1(i,j);
      u2t(i,j)=term1_2(i,j)+term2_2(i,j)+term3_2(i,j)...
              +term4_2(i,j)+term5_2(i,j);
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