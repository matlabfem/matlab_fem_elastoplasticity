function R=kelas3dstf(p,t,Young,nu)
%KELAS3DSTF Assembles the stiffness matrix for 3D linear elasticity 
% with P1 finite element  
%--------------------------------------------------------------------
% R=kelas3dstf(p,t,Young,nu)
%
% Input:
%      p   : Nodes coordinates, np*3
%      t   : Tetrahedron vertices, nt*4
%   Young  : Young modulus (scalar)
%      nu  : Poisson ratio (scalar)
%
% Output:
%      R   : Stiffness matrix, sparse (3*np)*(3*np)
%            symmetric, positive semidefinite
%--------------------------------------------------------------------
% (c) J. Koko, LIMOS 2015, koko@isima.fr
%-------------------------------------------------------------------- 
np=size(p,1); nn=3*np; nt=size(t,1);
mu=.5*Young/(1+nu); lam=Young*nu/((1+nu)*(1-2*nu));
C=zeros(6);C(1:3,1:3)=[lam+mu lam    lam;
                       lam    lam+mu lam;
                       lam    lam    lam+mu]; C=C+mu*eye(6);
% gradients of basis functions
[vol,g1,g2,g3,g4] = kpde3dgphi(p,t);   
B=cell(6,12); [B{:,:}]=deal(sparse(nt,1));
B(1,1:3:10)={g1(:,1) g2(:,1) g3(:,1) g4(:,1)};
B(2,2:3:11)={g1(:,2) g2(:,2) g3(:,2) g4(:,2)};
B(3,3:3:12)={g1(:,3) g2(:,3) g3(:,3) g4(:,3)};
B(4,1:3:10)={g1(:,2) g2(:,2) g3(:,2) g4(:,2)};
B(4,2:3:11)={g1(:,1) g2(:,1) g3(:,1) g4(:,1)};
B(5,1:3:10)={g1(:,3) g2(:,3) g3(:,3) g4(:,3)};
B(5,3:3:12)={g1(:,1) g2(:,1) g3(:,1) g4(:,1)};
B(6,2:3:11)={g1(:,3) g2(:,3) g3(:,3) g4(:,3)};
B(6,3:3:12)={g1(:,2) g2(:,2) g3(:,2) g4(:,2)};
E=cell(6,12);                          % E=C*B
for i=1:6
    for j=1:12
        E{i,j}=C(i,1)*B{1,j}+C(i,2)*B{2,j}+C(i,3)*B{3,j}...
              +C(i,4)*B{4,j}+C(i,5)*B{5,j}+C(i,6)*B{6,j};
    end
end
inodes=[1 1 1 2 2 2 3 3 3 4 4 4]; icomps=[1 2 3 1 2 3 1 2 3 1 2 3];
R=sparse(nn,nn); 
% under-diagonal entries
for i=1:12
    ik=3*(inodes(i)-1)+icomps(i); it=3*(t(:,inodes(i))-1)+icomps(i);
    for j=1:i-1
        jl=3*(inodes(j)-1)+icomps(j); jt=3*(t(:,inodes(j))-1)+icomps(j);
        Rij=B{1,ik}.*E{1,jl}+B{2,ik}.*E{2,jl}+B{3,ik}.*E{3,jl}...
           +B{4,ik}.*E{4,jl}+B{5,ik}.*E{5,jl}+B{6,ik}.*E{6,jl};
        R=R+sparse(it,jt,vol.*Rij,nn,nn);
    end
end
% transpose
R=R+R.';     
% transpose
for i=1:12                             % diagonal entries
    ik=3*(inodes(i)-1)+icomps(i); it=3*(t(:,inodes(i))-1)+icomps(i);
    Rij=B{1,ik}.*E{1,ik}+B{2,ik}.*E{2,ik}+B{3,ik}.*E{3,ik}...
       +B{4,ik}.*E{4,ik}+B{5,ik}.*E{5,ik}+B{6,ik}.*E{6,ik};
    R=R+sparse(it,it,vol.*Rij,nn,nn);
end
