function M=kpde3dmss(p,t,alfa)
%KPDE3DMSS Assembles the mass matrix with P1 finite element
%--------------------------------------------------------------------
% M=kpde3dmss(p,t,alfa)
%
%  Input:
%         p : Node coordinates, np*3
%         t : Tetrahedron vertices, nt*4  
%      alfa : scalar or column vector
%             If length(alfa)=np, alfa is defined on nodes
%             If length(alfa)=nt, alfa is defined on tetrahedra
%
%  Output:
%       M   : Mass matrix (symmetric positive definite) 
%             sparse np*np  
%--------------------------------------------------------------------
% (c) J. Koko, LIMOS 2015, koko@isima.fr
%--------------------------------------------------------------------
np=size(p,1); nt=size(t,1); [m1,m2]=size(alfa);
if (m1==1 || m1==nt)
    c=alfa; 
else
    c=(alfa(t(:,1))+alfa(t(:,2))+alfa(t(:,3))+alfa(t(:,4)))/4;
end
% elements area
vol=kpde3dgphi(p,t); c=c.*vol/20;  
M=sparse(np,np);
% under-diagonal entries
for i=1:4
    for j=1:i-1             
        M=M+sparse(t(:,i),t(:,j),c,np,np);
    end
end
% transpose
M=M+M.';                           
% diagonal entries
for i=1:4
    M=M+sparse(t(:,i),t(:,i),2*c,np,np); 
end
