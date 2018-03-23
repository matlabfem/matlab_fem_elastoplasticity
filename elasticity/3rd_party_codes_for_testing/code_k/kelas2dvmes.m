function S=kelas2dvmes(p,t,u,E,nu)
%KELAS2DVMES Computes 2D Von Mises effective stress
%--------------------------------------------------------------------
%      S=kelas2dvmes(p,t,u,E,nu)           
%
% Input:
%         p : node coodinates, np∗2
%         t : triangle vertices, nt∗3
%         u : displacements field, (2*np)*1
%      E,nu : Young modulus and Poisson ratio
%      S    : Von Mises effective stress, np*1
%--------------------------------------------------------------------
% (c) J. Koko, LIMOS 2006-2015, koko@isima.fr
%--------------------------------------------------------------------
np=size(p,1); nn=2*np; lam=E*nu/((1+nu)*(1-2*nu));  mu=E/(2*(1+nu));
it1=t(:,1); it2=t(:,2); it3=t(:,3);
%  elements area and gradient of basis functions
[ar,phi1,phi2,phi3]=kpde2dgphi(p,t);
% strain
u1=u(1:2:end); u2=u(2:2:end); 
uh=[u1(it1) u1(it2) u1(it3)]; vh=[u2(it1) u2(it2) u2(it3)];
e11=uh(:,1).*phi1(:,1)+uh(:,2).*phi2(:,1)+uh(:,3).*phi3(:,1);
e22=vh(:,1).*phi1(:,2)+vh(:,2).*phi2(:,2)+vh(:,3).*phi3(:,2);
e12=uh(:,1).*phi1(:,2)+uh(:,2).*phi2(:,2)+uh(:,3).*phi3(:,2)...
   +vh(:,1).*phi1(:,1)+vh(:,2).*phi2(:,1)+vh(:,3).*phi3(:,1);
sig11=(lam+2*mu)*e11+lam*e22; sig22=lam*e11+(lam+2*mu)*e22; sig12=mu*e12;
clear uh vh e11 e22 e12 
% patch
arp=sparse(it1,1,ar,np,1)+sparse(it2,1,ar,np,1)+sparse(it3,1,ar,np,1);
sm1=ar.*sig11; sm2=ar.*sig22; sm12=ar.*sig12;
s1=sparse(it1,1,sm1,np,1)+sparse(it2,1,sm1,np,1)+sparse(it3,1,sm1,np,1); 
s2=sparse(it1,1,sm2,np,1)+sparse(it2,1,sm2,np,1)+sparse(it3,1,sm2,np,1); 
s12=sparse(it1,1,sm12,np,1)+sparse(it2,1,sm12,np,1)+sparse(it3,1,sm12,np,1); 
s1=full(s1./arp); s2=full(s2./arp); s12=full(s12./arp);
% Von mises effective stress
delta=sqrt((s1-s2).^2+4*s12.^2); sp1=s1+s2+delta; sp2=s1+s2-delta;
S=sqrt(sp1.^2+sp2.^2-sp1.*sp2);