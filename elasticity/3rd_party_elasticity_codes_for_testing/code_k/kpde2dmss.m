function M=kpde2dmss(p,t,alfa)
%KPDE2DMSS Assembles the mass matrix with P1 finite element
%--------------------------------------------------------------------
% M=kpde2dmss(p,t,alfa)
%
%  Input:
%         p : Node coordinates, np*2
%         t : Triangle vertices, nt*3  
%      alfa : scalar or column vector
%             If length(alfa)=np, alfa is defined on nodes
%             If length(alfa)=nt, alfa is defined on triangles
%
%  Output:
%       M   : mass  matrix (symmetric positive definite) 
%             sparse np*np  
%--------------------------------------------------------------------
% (c) J. Koko, LIMOS 2006-2015, koko@isima.fr
%--------------------------------------------------------------------
np=size(p,1); nt=size(t,1); [m1,m2]=size(alfa);
if (m1==1 || m1==nt)
    c=alfa; 
else
    c=(alfa(t(:,1))+alfa(t(:,2))+alfa(t(:,3)))/3; 
end
ar=kpde2dgphi(p,t); c=c.*ar/12;   
M=sparse(np,np);
% under-diagonal entries
for i=1:3
    for j=1:i-1           
        M=M+sparse(t(:,i),t(:,j),c,np,np);
    end
end
M=M+M.';  
% diagonal entries
for i=1:3
    M=M+sparse(t(:,i),t(:,i),2*c,np,np); 
end
