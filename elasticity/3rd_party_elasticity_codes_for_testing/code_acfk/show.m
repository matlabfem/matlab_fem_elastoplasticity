function show(elements,dirichlet,neumann,coordinates,u,lambda,mu)
%SHOW  Plots three-dimensional solution
%    SHOW(ELEMENTS,DIRICHLET,NEUMANN,COORDINATES,U,LAMBDA,MU) plots the
%    strained mesh and visualizes the stresses in grey tones.
%
%
%   See also FEM_LAME3D.

%    J. Alberty, C. Carstensen, S. A. Funken, and R. Klose  07-03-00
%    File <show.m> in $(HOME)/acfk/fem_lame3d/angle/

E=zeros(4*size(elements,1),3);
 C=zeros(size(elements,1),1);
 for j=1:size(elements,1)
   PhiGrad=[1,1,1,1;coordinates(elements(j,:),:)']\[zeros(1,3);eye(3)]; 
   U_Grad = u([1;1;1]*3*elements(j,:)-[2;1;0]*[1,1,1,1])*PhiGrad;
   SIGMA = lambda * trace(U_Grad)*eye(3)+mu*(U_Grad + U_Grad');
   C(j) = sqrt(sum(eig(SIGMA).^2));   
 end
 Area = zeros(size(elements,1),1);
 AreaOmega = zeros(max(max(elements)),1);
 AvC = zeros(max(max(elements)),1);
 for j=1:size(elements,1)
   Area = det([1,1,1,1;coordinates(elements(j,:),:)'])/6;
   AreaOmega(elements(j,:))=AreaOmega(elements(j,:))+Area;
   AvC(elements(j,:)) = AvC(elements(j,:))+Area*[1;1;1;1]*C(j);
 end
 AvC = AvC./AreaOmega;
 E=[dirichlet;neumann];
  factor=100.0;
  %colormap(1-gray);
  trisurf(E,factor*u(1:3:size(u,1))+coordinates(:,1), ...
        factor*u(2:3:size(u,1))+coordinates(:,2), ...
        factor*u(3:3:size(u,1))+coordinates(:,3),AvC,'facecolor','interp');
  view(-50,35)
  axis equal; axis([-10 70 -100 20 0 40])
  xlabel('x'); ylabel('y'); zlabel('z');
  
