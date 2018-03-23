function b=kelas3drhsn(p,fneum,g1,g2,g3)
%KPDE3DRHSN Assembles the RHS contribution of Neumann boundary conditions.
%--------------------------------------------------------------------
% b=kelas3drhsn(p,fneum,g1,g2,g3)
%
% input:
%        p : Nodes coordinates, np*3
%    fneum : Boundary faces (triangles), nf*3  
% g1,g2,g3 : Neumann boundary condition,column vectors np*1
%            g ~=0 on Neumann nodes 
%
% Output:
%        b : Right-hand side, np*1  
%--------------------------------------------------------------------
% (c) J. Koko, LIMOS 2015, koko@isima.fr
%--------------------------------------------------------------------
np=size(p,1); nt=size(fneum,1);

% faces area
it1=fneum(:,1); it2=fneum(:,2); it3=fneum(:,3); 
v1=p(it2,:)-p(it1,:); v2=p(it3,:)-p(it1,:);
v=cross(v1,v2,2);  ar=sqrt(sum(v.*v,2))/2;

% g at centers of mass
gh1=ar.*(g1(it1)+g1(it2)+g1(it3))/9; 
gh2=ar.*(g2(it1)+g2(it2)+g2(it3))/9; 
gh3=ar.*(g3(it1)+g3(it2)+g3(it3))/9; 
 
% assembly
gg1=sparse(it1,1,gh1,np,1)+sparse(it2,1,gh1,np,1)+sparse(it3,1,gh1,np,1);
gg2=sparse(it1,1,gh2,np,1)+sparse(it2,1,gh2,np,1)+sparse(it3,1,gh2,np,1);
gg3=sparse(it1,1,gh3,np,1)+sparse(it2,1,gh3,np,1)+sparse(it3,1,gh3,np,1);
b=zeros(3*np,1);
b(1:3:3*np-2)=full(gg1); b(2:3:3*np-1)=full(gg2); b(3:3:3*np)=full(gg3);