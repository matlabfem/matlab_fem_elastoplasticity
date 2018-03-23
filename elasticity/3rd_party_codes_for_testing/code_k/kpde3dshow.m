function []=kpde3dshow(p,t,u,it)
%KPDE3DSHOW  KPDE3D Visualization function
%--------------------------------------------------------------------
%       kpde3dshow(p,t,u)
%       kpde3dshow(p,t,u,it)
%
% Input:
%      p : node coodinate, np∗3
%      t : tetrahedron vertices, nt∗4
%      u : PDE solution,  np*1 
%     it : set of tetrahedra defining a subdomain of interest
%--------------------------------------------------------------------
% (c) J. Koko, LIMOS 2015, koko@isima.fr
%--------------------------------------------------------------------
ii=[1:size(t,1)]';
if (nargin == 4)
    ii=it;
end
% faces of tetrahedra
fh0=[t(ii,[1,2,3]);
   t(ii,[1,2,4]);
   t(ii,[1,3,4]);
   t(ii,[2,3,4])];
fh=sort(fh0,2);
[~,ie,je]=unique(fh,'rows');
% boundary faces 
f=fh(ie,:); 
nf=size(f,1);
hf=accumarray(je,1,[nf 1]); ib=find(hf==1);
tb=f(ib,:);
% visualization with trisurf
trisurf(tb,p(:,1),p(:,2),p(:,3),u,'facecolor','interp'), axis image