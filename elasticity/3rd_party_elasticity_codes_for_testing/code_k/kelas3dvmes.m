function S=kelas3dvmes(p,t,u,E,nu)
%KELAS3DVMES Computes 3D Von Mises effective stress 
%--------------------------------------------------------------------
%      S=kelas3dvmes(p,t,u,E,nu)   
%
% Input:
%         p : node coodinate, np∗3
%         t : tetrahedron vertices, nt∗4
%         u : displacements field, (3*np)*1
%      E,nu : Young modulus and Poisson ratio
%      S    : Von Mises effective stress, np*1
%--------------------------------------------------------------------
% (c) J. Koko, LIMOS 2015, koko@isima.fr
%--------------------------------------------------------------------
np=size(p,1); nn=3*np; lam=E*nu/((1+nu)*(1-2*nu));  mu=E/(2*(1+nu));
C=zeros(6); C(1:3,1:3)=[lam+mu lam lam; lam lam+mu lam;lam lam lam+mu];
C=C+mu*eye(6);
% elements volume and gradient of basis functions 
[vol,g1,g2,g3,g4]=kpde3dgphi(p,t);
% deformations 
it1=t(:,1); it2=t(:,2); it3=t(:,3); it4=t(:,4);
u1=u(1:3:nn-2); u2=u(2:3:nn-1);  u3=u(3:3:nn);
e11=u1(it1).*g1(:,1)+u1(it2).*g2(:,1)+u1(it3).*g3(:,1)+u1(it4).*g4(:,1);
e22=u2(it1).*g1(:,2)+u2(it2).*g2(:,2)+u2(it3).*g3(:,2)+u2(it4).*g4(:,2);
e33=u3(it1).*g1(:,3)+u3(it2).*g2(:,3)+u3(it3).*g3(:,3)+u3(it4).*g4(:,3);
e12=u1(it1).*g1(:,2)+u1(it2).*g2(:,2)+u1(it3).*g3(:,2)+u1(it4).*g4(:,2)...
   +u2(it1).*g1(:,1)+u2(it2).*g2(:,1)+u2(it3).*g3(:,1)+u2(it4).*g4(:,1);
e13=u1(it1).*g1(:,3)+u1(it2).*g2(:,3)+u1(it3).*g3(:,3)+u1(it4).*g4(:,3)...
   +u3(it1).*g1(:,1)+u3(it2).*g2(:,1)+u3(it3).*g3(:,1)+u3(it4).*g4(:,1);
e23=u2(it1).*g1(:,3)+u2(it2).*g2(:,3)+u2(it3).*g3(:,3)+u2(it4).*g4(:,3)...
   +u3(it1).*g1(:,2)+u3(it2).*g2(:,2)+u3(it3).*g3(:,2)+u3(it4).*g4(:,2);
% stress
s=cell(6,1); s{4}=mu*vol.*e12; s{5}=mu*vol.*e13; s{6}=mu*vol.*e23;
for i=1:3
    s{i}=vol.*(C(i,1)*e11+C(i,2)*e22+C(i,3)*e33);
end
clear e11 e22 e33 e13 e23 e12
% patch
volp=sparse(it1,1,vol,np,1)+sparse(it2,1,vol,np,1)...
    +sparse(it3,1,vol,np,1)+sparse(it4,1,vol,np,1);
% average stress
for i=1:6
    sm{i}=sparse(it1,1,s{i},np,1)+sparse(it2,1,s{i},np,1)...
         +sparse(it3,1,s{i},np,1)+sparse(it4,1,s{i},np,1);
    sm{i}=full(sm{i}./volp);
end
% Von Mises effective stress
I1=sm{1}+sm{2}+sm{3};
I2=sm{1}.*sm{2}+sm{2}.*sm{3}+sm{3}.*sm{1}-sm{4}.^2-sm{5}.^2-sm{6}.^2;
S=full(sqrt(I1.^2-3*I2));

