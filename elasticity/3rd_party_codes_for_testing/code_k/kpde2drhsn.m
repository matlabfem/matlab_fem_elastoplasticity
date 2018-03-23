function b=kpde2drhsn(p,eneum,g)
%KPDE2DRHSN Assembles the RHS contribution of Neumann boundary conditions.
%--------------------------------------------------------------------
% b=edp2drhsn(p,eneum,g)
%
% input:
%        p : Nodes coordinates, np*2
%    eneum : Boundary edges, ne*2  
%        g : Neumann boundary condition, np*1
%            g ~=0 on Neumann nodes 
%
% Output:
%        b : Right-hand side, np*1  
%--------------------------------------------------------------------
% (c) J. Koko, LIMOS 2006-2015, koko@isima.fr
%--------------------------------------------------------------------
np=size(p,1); ie1=eneum(:,1); ie2=eneum(:,2);
% segments length
xy=p(ie2,:)-p(ie1,:); ee=sqrt(sum(xy.^2,2));  
gg=ee.*(g(ie1)+g(ie2))/2;    
b=sparse(ie1,1,gg,np,1)+sparse(ie2,1,gg,np,1); 
b=full(b);
