function draw_mesh(coord,surf,elem_type)

% =========================================================================
%
%  This function draws mesh and nodal point on the surface of the body
%
%  input data:
%    coord - coordinates of the nodes, size(coord)=(3,n_n) where n_n is a
%            number of nodes
%    surf  - array containing numbers of nodes defining each surface element,
%            size(surf)=(n_p,n_s), n_s = number of surface elements
%    elem_type - the type of finite elements; available choices:
%                'P1', 'P2', 'Q1', 'Q2'
%
% ======================================================================
%

  figure
  hold on
  if strcmp(elem_type,'P1')||strcmp(elem_type,'P2')
      patch('Faces',surf(1:3,:)','Vertices',coord','FaceVertexCData',...
           0*ones(size(coord,2),1),'FaceColor','white','EdgeColor','blue'); 
  else
      patch('Faces',surf(1:4,:)','Vertices',coord','FaceVertexCData',...
           0*ones(size(coord,2),1),'FaceColor','white','EdgeColor','blue');
  end
  ind=unique(surf(:));
  plot3( coord(1,ind),coord(2,ind),coord(3,ind), 'b.', 'MarkerSize',10);
  axis equal;  % real ratios
  view(3);     % standard view ve 3D
  hold off;
  axis off;
end