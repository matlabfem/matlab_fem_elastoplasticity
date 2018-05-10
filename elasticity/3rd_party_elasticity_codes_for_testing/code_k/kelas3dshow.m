function kelas3dshow(p,t,u,mag,sf)
%KELAS3DSHOW 3D linear elasticity visualization
%--------------------------------------------------------------------
% kelas3dshow(p,t,u,mag)       - deformed mesh
% kelas3dshow(p,t,u,mag,sf)    - deformed mesh + sf 
%
% Input:
%      p  : node coordinates, np*3
%      t  : tetrahedron vertices, nt*4
%      u  : Displacements vector, (3*np)*1
%    mag  : Magnification factor, positive sclar 
%      s  : Variable to be visualized (e.g. shear energy density), np*1
%
%--------------------------------------------------------------------
% (c) J. Koko, LIMOS 2015, koko@isima.fr
%-------------------------------------------------------------------- 
np=size(p,1);  pu=zeros(np,3); 
pu(:,1)=p(:,1)+mag*u(1:3:end); 
pu(:,2)=p(:,2)+mag*u(2:3:end);
pu(:,3)=p(:,3)+mag*u(3:3:end);
 
if (nargin==4) 
    kpde3dshow(pu,t,zeros(np,1));  % deformed mesh
    axis on
else
   kpde3dshow(pu,t,sf);            % deformed mesh + stress
   axis on, colorbar('South')
end