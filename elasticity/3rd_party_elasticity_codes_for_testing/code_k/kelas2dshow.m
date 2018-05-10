function kelas2dshow(p,t,u,mag,sf)
%KELAS2DSHOW Visualisation du maillage déformé et/ou des contraintes
 %--------------------------------------------------------------------
np=size(p,1);  pu=zeros(np,2); 
pu(:,1)=p(:,1)+mag*u(1:2:end); pu(:,2)=p(:,2)+mag*u(2:2:end);
 
% 2. maillage déformé + contraintes en niveau de gris 
trisurf(t,pu(:,1),pu(:,2),zeros(np,1),sf,'facecolor','interp')
colorbar('East'), view(2), return
 
 