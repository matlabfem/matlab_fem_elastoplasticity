function R=kpde2dstf(p,t,nu)
%KPDE2DSTF Assembles the stiffness matrix with P1 finite element
%--------------------------------------------------------------------
% R=kpde2dstf(p,t,nu)
%
%  Input:
%         p : Node coordinates, np*2
%         t : Triangle vertices, nt*3  
%        nu : scalar or column vector
%             If length(nu)=np, nu is defined on nodes
%             If length(nu)=nt, nu is defined on triangles
%
%  Output:
%       R  : Stiffness matrix (symmetric positive semidefinite) 
%            sparse np*np  
%--------------------------------------------------------------------
% (c) J. Koko, LIMOS 2006-2015, koko@isima.fr
%--------------------------------------------------------------------
np=size(p,1); nt=size(t,1); [m1,m2]=size(nu);
if (m1==1 || m1==nt)
    c=nu; 
else
    c=(nu(t(:,1))+nu(t(:,2))+nu(t(:,3)))/3; 
end
[ar,g1,g2,g3]=kpde2dgphi(p,t);  c=c.*ar;
% cell-array of gradients    
g={g1 g2 g3};       
R=sparse(np,np);
% under-diagonal entries
for i=1:3
    for j=1:i-1       
        R=R+sparse(t(:,i),t(:,j),c.*sum(g{i}.*g{j},2),np,np);
    end
end
R=R+R.';                     
% diagonal entries
for i=1:3
    R=R+sparse(t(:,i),t(:,i),c.*sum(g{i}.*g{i},2),np,np);  
end
