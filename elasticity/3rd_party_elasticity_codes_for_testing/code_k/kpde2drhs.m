function b=kpde2drhs(p,t,f)
%KPDE2DRHS Assembles the right-hand side with P1 finite element
%--------------------------------------------------------------------
% b=kpde2drhs(p,t,f)
%
% Input:
%      p : Nodes coordinates, np*2
%      t : Triangle vertices, nt*3  
%      f : Source term, column vector np*1 or nt*1 
%          If length(f)=np, f is defined on nodes
%          If length(f)=nt, f is defined on triangles
%
% Output:
%     b : Right-hand side, np*1  
%--------------------------------------------------------------------
% (c) J. Koko, LIMOS 2006-2015, koko@isima.fr
%--------------------------------------------------------------------
np=size(p,1); nt=size(t,1); [m1,m2]=size(f);
if (m1==1 || m1==nt)
    ff=f;
else
    ff=(f(t(:,1))+f(t(:,2))+f(t(:,3)))/3; 
end
ar=kpde2dgphi(p,t); ff=ff.*ar/3;
b=sparse(np,1); 
for i=1:3
    b=b+sparse(t(:,i),1,ff,np,1); 
end
b=full(b);
