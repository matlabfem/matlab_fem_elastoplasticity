function [vol,g1,g2,g3,g4]=kpde3dgphi(p,t)
%KPDE3DGPHI Elenemts volume and gradient of basis functions
%--------------------------------------------------------------------
%            vol=kpde3dgphi(p,t) 
% [vol,g1,g2,g3]=kpde3dgphi(p,t)
%
%  Input:
%            p : Node coordinates, np*3
%            t : Tetrahedron vertices, nt*4  
%
%  Output:
%          vol : elements volume, nt*1
%  g1,g2,g3,g4 : gradient of basis functions, nt*3
%--------------------------------------------------------------------
% (c) J. Koko, LIMOS 2015, koko@isima.fr
%--------------------------------------------------------------------
it1=t(:,1); it2=t(:,2); it3=t(:,3); it4=t(:,4);
x1=p(it1,1); x2=p(it2,1); x3=p(it3,1); x4=p(it4,1);
y1=p(it1,2); y2=p(it2,2); y3=p(it3,2); y4=p(it4,2);
z1=p(it1,3); z2=p(it2,3); z3=p(it3,3); z4=p(it4,3);
% 3x3 Jacobian matrix & determinant
J11=x2-x1; J12=y2-y1; J13=z2-z1; 
J21=x3-x1; J22=y3-y1; J23=z3-z1;
J31=x4-x1; J32=y4-y1; J33=z4-z1;
det=J11.*(J22.*J33-J32.*J23)+J12.*(J31.*J23-J21.*J33)+J13.*(J21.*J32-J31.*J22);
% elements volume
vol=abs(det)/6;                     
if (nargout == 1), return; end
% Jacobian inverse C
C11=(J22.*J33-J32.*J23)./det; C12=(J13.*J32-J12.*J33)./det; C13=(J12.*J23-J13.*J22)./det;
C21=(J31.*J23-J21.*J33)./det; C22=(J11.*J33-J13.*J31)./det; C23=(J21.*J13-J23.*J11)./det;
C31=(J21.*J32-J31.*J22)./det; C32=(J12.*J31-J32.*J11)./det; C33=(J11.*J22-J12.*J21)./det;
% gradients of basis functions
g2=[C11 C21 C31]; g3=[C12 C22 C32]; g4=[C13 C23 C33]; g1=-g2-g3-g4;
