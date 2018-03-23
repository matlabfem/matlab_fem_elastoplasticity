function b=kelas2drhsn(p,eneum,g1,g2)
%KELAS2DRHSN Assembles the RHS contribution of Neumann boundary conditions.
%--------------------------------------------------------------------
% b=kelas2drhsn(p,eneum,gx,gy)
%
% input:
%        p : Nodes coordinates, np*2
%    fneum : Boundary edges, ne*2  
%    g1,g2 : Neumann boundary condition, np*1
%            g ~=0 on Neumann nodes 
%
% Output:
%        b : Right-hand side, 2*np*1  
%--------------------------------------------------------------------
% (c) J. Koko, LIMOS 2006-2015, koko@isima.fr
%--------------------------------------------------------------------
np=size(p,1); ni=size(eneum,1);

% edges length
xy=p(eneum(:,2),:)-p(eneum(:,1),:); le=sqrt(sum(xy.^2,2));

% g at centers of mass
gh1=(g1(eneum(:,1))+g1(eneum(:,2))).*le/4; 
gh2=(g2(eneum(:,1))+g2(eneum(:,2))).*le/4; 
 
% assembly
b=zeros(2*np,1);
b(1:2:end)=full(sparse(eneum(:,1),1,gh1,np,1)+sparse(eneum(:,2),1,gh1,np,1));
b(2:2:end)=full(sparse(eneum(:,1),1,gh2,np,1)+sparse(eneum(:,2),1,gh2,np,1));