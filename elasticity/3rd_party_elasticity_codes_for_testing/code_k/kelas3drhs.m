function b=kelas3drhs(p,t,fx,fy,fz)
%KELAS3DRHS Assembles the right-hand side with P1 finite element
%--------------------------------------------------------------------
% b=kpde3drhs(p,t,fx,fy,fz)
%
% Input:
%      p : Nodes coordinates, np*3
%      t : Tetrahedron vertices, nt*4  
%      f : Source term, column vector np*1  
%
% Output:
%      b : Right-hand side, (3*np)*1  
%--------------------------------------------------------------------
% (c) J. Koko, LIMOS 2015, koko@isima.fr
%--------------------------------------------------------------------
np=size(p,1); nt=size(t,1); nn=3*np;
it1=t(:,1); it2=t(:,2); it3=t(:,3); it4=t(:,4);

% f at centers of mass
f1=(fx(it1)+fx(it2)+fx(it3)+fx(it4))/4; 
f2=(fy(it1)+fy(it2)+fy(it3)+fy(it4))/4; 
f3=(fz(it1)+fz(it2)+fz(it3)+fz(it4))/4; 
 
% elements volume
vol= kpde3dgphi(p,t);   
 
% assembly
f=[f1.*vol/4 f2.*vol/4 f3.*vol/4]; 
ff=sparse(np,3);
for i=1:4
    ff(:,1)=ff(:,1)+sparse(t(:,i),1,f(:,1),np,1);
    ff(:,2)=ff(:,2)+sparse(t(:,i),1,f(:,2),np,1);
    ff(:,3)=ff(:,3)+sparse(t(:,i),1,f(:,3),np,1);
end
ff=full(ff); b=zeros(nn,1);
b(1:3:nn-2)=ff(:,1); b(2:3:nn-1)=ff(:,2); b(3:3:nn)=ff(:,3);
