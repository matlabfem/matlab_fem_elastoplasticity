function [ar,g1,g2,g3]=kpde2dgphi(p,t)
%KPDE2DGPHI Elements area and gradient of basis functions
%--------------------------------------------------------------------
%            ar=kpde2dgphi(p,t) 
% [ar,g1,g2,g3]=kpde2dgphi(p,t)
%
%  Input:
%         p : Node coordinates, np*2
%         t : Triangle vertices, nt*3  
%
%  Output:
%        ar : elements area, nt*1
%  g1,g2,g3 : gradient of basis functions, nt*2
%--------------------------------------------------------------------
% (c) J. Koko, LIMOS 2006-2015, koko@isima.fr
%--------------------------------------------------------------------
np=size(p,1); nt=size(t,1);
x21=p(t(:,2),1)-p(t(:,1),1); y21=p(t(:,2),2)-p(t(:,1),2); 
x32=p(t(:,3),1)-p(t(:,2),1); y32=p(t(:,3),2)-p(t(:,2),2);
x31=p(t(:,3),1)-p(t(:,1),1); y31=p(t(:,3),2)-p(t(:,1),2);
% triangles area
ar=(x21.*y31-y21.*x31)/2;  
if (nargout==1), return; end
% gradients of basis functions
g1=.5*[-y32./ar x32./ar]; g2=.5*[y31./ar -x31./ar]; g3=.5*[-y21./ar x21./ar]; 
