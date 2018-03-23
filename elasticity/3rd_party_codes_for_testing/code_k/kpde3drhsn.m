function gh=kpde3drhsn(p,fneum,g)
%KPDE3DRHSN Assembles the RHS contribution of Neumann boundary conditions.
%--------------------------------------------------------------------
% b=edp3drhsn(p,fneum,g)
%
% input:
%        p : Nodes coordinates, np*3
%    fneum : Boundary faces (triangles), nf*3  
%        g : Neumann boundary condition, np*1
%            g ~=0 on Neumann nodes 
%
% Output:
%        b : Right-hand side, np*1  
%--------------------------------------------------------------------
% (c) J. Koko, LIMOS 2015, koko@isima.fr
%--------------------------------------------------------------------
%--------------------------------------------------------------------
np=size(p,1); it1=fneum(:,1); it2=fneum(:,2); it3=fneum(:,3);
% faces area
v1=p(it2,:)-p(it1,:); v2=p(it3,:)-p(it1,:); v=cros(v1,v2,2);
farea=sqrt(sum(v.*v,2))/2;                 
gg=farea.*(g(it1)+g(it2)+g(it3))/3; 
b=sparse(it1,1,gg,np,1)+sparse(it2,1,gg,np,1)+sparse(it3,1,gg,np,1);
b=full(b);
