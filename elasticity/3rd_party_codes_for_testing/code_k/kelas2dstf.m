function R=kelas2dstf(p,t,Young,nu)
%KELAS2DSTF Assembles the stiffness matrix for 2D linear elasticity 
% with P1 finite element  
%--------------------------------------------------------------------
% R=kelas2dstf(p,t,Young,nu)
%
% Input:
%      p   : Nodes coordinates, np*2
%      t   : Triangle vertices, nt*3
%   Young  : Young modulus (scalar)
%      nu  : Poisson ratio (scalar)
%
% Output:
%      R   : Stiffness matrix, sparse (2*np)*(2*np)
%            symmetric, positive semidefinite
%--------------------------------------------------------------------
% (c) J. Koko, LIMOS 2015, koko@isima.fr
%-------------------------------------------------------------------- 
np=size(p,1); nn=2*np; nt=size(t,1);
mu=.5*Young/(1+nu); lam=Young*nu/((1+nu)*(1-2*nu));
C=[lam+2*mu lam       0; 
   lam      lam+2*mu  0; 
   0        0         mu];
% gradients of basis functions
[ar,g1,g2,g3]=kpde2dgphi(p,t);         
B=cell(3,6); [B{:,:}]=deal(sparse(nt,1));    
B{1,1}=g1(:,1); B{1,3}=g2(:,1); B{1,5}=g3(:,1);
B{2,2}=g1(:,2); B{2,4}=g2(:,2); B{2,6}=g3(:,2);
B{3,1}=g1(:,2); B{3,2}=g1(:,1); B{3,3}=g2(:,2); 
B{3,4}=g2(:,1); B{3,5}=g3(:,2); B{3,6}=g3(:,1);
E=cell(3,6);                           
for i=1:3
    for j=1:6                    
         E{i,j}=C(i,1)*B{1,j}+C(i,2)*B{2,j}+C(i,3)*B{3,j};
    end
end
inodes=[1 1 2 2 3 3]; icomps=[1 2 1 2 1 2];
R=sparse(nn,nn);                      
% under-diagonal entries
for i=1:6
    ik=2*(inodes(i)-1)+icomps(i); it=2*(t(:,inodes(i))-1)+icomps(i);
    for j=1:i-1
        jl=2*(inodes(j)-1)+icomps(j); jt=2*(t(:,inodes(j))-1)+icomps(j);
        Rij=B{1,ik}.*E{1,jl}+B{2,ik}.*E{2,jl}+B{3,ik}.*E{3,jl};
        R=R+sparse(it,jt,ar.*Rij,nn,nn);
    end
end
% transpose
R=R+R.';      
% diagonal entries
for i=1:6                              
    ik=2*(inodes(i)-1)+icomps(i); it=2*(t(:,inodes(i))-1)+icomps(i);
    Rij=B{1,ik}.*E{1,ik}+B{2,ik}.*E{2,ik}+B{3,ik}.*E{3,ik};
    R=R+sparse(it,it,ar.*Rij,nn,nn);
end
