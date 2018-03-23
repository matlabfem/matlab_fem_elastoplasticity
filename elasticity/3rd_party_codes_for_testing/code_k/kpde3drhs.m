function b=kpde3drhs(p,t,f)
%KPDE3DRHS Assembles the right-hand side with P1 finite element
%--------------------------------------------------------------------
% b=kpde3drhs(p,t,f)
%
% Input:
%      p : Nodes coordinates, np*3
%      t : Tetrahedron vertices, nt*4  
%      f : Source term, column vector np*1 or nt*1 
%          If length(f)=np, f is defined on nodes
%          If length(f)=nt, f is defined on tetrahedra
%
% Output:
%      b : Right-hand side, np*1  
%--------------------------------------------------------------------
% (c) J. Koko, LIMOS 2015, koko@isima.fr
%--------------------------------------------------------------------
 
np=size(p,1); nt=size(t,1); [m1,m2]=size(f);
if ( m1==1 || m1==nt ) 
    ff=f ;
else
    ff=(f(t(:,1))+f(t(:,2))+f(t(:,3))+f(t(:,4)))/4; 
end
vol=kpde3dgphi(p,t);  ff=ff.*vol/4;
b=sparse(np,1); 
for i=1:4
    b=b+sparse(t(:,i),1,ff,np,1); 
end
b=full(b);
