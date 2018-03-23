function R=kpde3dstf(p,t,nu)
%KPDE3DSTF Assembles the stiffness matrix with P1 finite element
%--------------------------------------------------------------------
% R=kpde3dstf(p,t,nu)
%
%  Input:
%         p : Node coordinates, np*3
%         t : Tetrahedron vertices, nt*4  
%        nu : scalar or column vector
%             If length(nu)=np, nu is defined on nodes
%             If length(nu)=nt, nu is defined on tetrahedra
%
%  Output:
%       R  : Stiffness matrix (symmetric positive semidefinite) 
%            sparse np*np  
%--------------------------------------------------------------------
% (c) J. Koko, LIMOS 2015, koko@isima.fr
%--------------------------------------------------------------------
np=size(p,1); nt=size(t,1); [m1,m2]= size(nu) ;
if (m1==1 || m1==nt)
    c=nu ; 
else
    c=(nu(t(:,1))+nu(t(:,2))+nu(t(:,3))+nu(t(:,4)))/4; 
end
% gradients of basis functions
[vol,g1,g2,g3,g4] = kpde3dgphi(p,t);   
c=c.*vol; 
% cell-array of gradient
g={g1 g2 g3 g4};             
R=sparse(np,np);
% under-diagonal entries
for i=1:4
    for j=1:i-1                
        R=R+sparse(t(:,i),t(:,j),c.*sum(g{i}.*g{j},2),np,np);
    end
end
% transpose
R=R+R.';   
% diagonal entries
for i=1:4
    R=R+sparse(t(:,i),t(:,i),c.*sum(g{i}.^2,2),np,np); 
end
