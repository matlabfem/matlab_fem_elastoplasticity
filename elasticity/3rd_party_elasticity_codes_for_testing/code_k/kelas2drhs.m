function b=kelas2drhs(p,t,fx,fy)
%KELAS2DRHS Assembles the right-hand side with P1 finite element
%--------------------------------------------------------------------
% b=kelas2drhs(p,t,fx,fy) 
%
% Input:
%      p : Nodes coordinates, np*2
%      t : Triangle vertices, nt*3  
%      f : Source term, column vector np*1  
%
% Output:
%      b : Right-hand side, np*1  
%-------------------------------------------------------------------- 
% (c) J. Koko, LIMOS 2006-2015, koko@isima.fr
%--------------------------------------------------------------------
np=size(p,1); nt=size(t,1); nn=2*np;

% f at centers of mass
f1=(fx(t(:,1))+fx(t(:,2))+fx(t(:,3)))/3; 
f2=(fy(t(:,1))+fy(t(:,2))+fy(t(:,3)))/3; 

% triangles area
area=kpde2dgphi(p,t);

% assembly
f1=f1.*area/3; f2=f2.*area/3;
b=zeros(nn,1); 
b(1:2:nn)=full(sparse(t(:,1),1,f1,np,1)+sparse(t(:,2),1,f1,np,1)...
          +sparse(t(:,3),1,f1,np,1));
b(2:2:nn)=full(sparse(t(:,1),1,f2,np,1)+sparse(t(:,2),1,f2,np,1)...
          +sparse(t(:,3),1,f2,np,1)); 