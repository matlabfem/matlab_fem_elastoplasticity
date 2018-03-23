function draw_mesh(coord,elem,elem_type)

% =========================================================================
%
%  This function draws mesh and nodal point on the surface of the body
%
%  input data:
%    coord - coordinates of the nodes, size(coord)=(2,n_n) where n_n is a
%            number of nodes
%    surf - n_p x n_s array containing numbers of nodes defining each
%           surface element, n_s = number of surface elements
%    elem_type - the type of finite elements; available choices:
%                'P1', 'P2', 'Q1', 'Q2'
%
% ======================================================================
%

  figure
  hold on
  coord_aux = [coord;zeros(1,size(coord,2))];
  if strcmp(elem_type,'P1')||strcmp(elem_type,'P2')      
    patch('Faces',elem(1:3,:)','Vertices',coord_aux','FaceVertexCData',...
           0*ones(size(coord,2),1),'FaceColor','white','EdgeColor','blue'); 
  else
    patch('Faces',elem(1:4,:)','Vertices',coord_aux','FaceVertexCData',...
           0*ones(size(coord,2),1),'FaceColor','white','EdgeColor','blue');
  end
%   ind=unique(surf(:));
  plot( coord(1,:),coord(2,:), 'b.', 'MarkerSize',10);
  axis equal;  %realne pomery
  view(2);     %standartni pohled ve 2D
  hold off;
  axis off;
end