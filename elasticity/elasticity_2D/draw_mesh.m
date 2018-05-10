function draw_mesh(coord,elem,elem_type)

% =========================================================================
%
%  This function draws mesh and nodal point on the surface of the body
%
%  input data:
%    coord - coordinates of the nodes, size(coord)=(2,n_n) where n_n is a
%            number of nodes
%    elem  - array containing numbers of nodes defining each element, 
%            size(elem)=(n_p, n_e), n_e = number of elements
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

  plot( coord(1,:),coord(2,:), 'b.', 'MarkerSize',10);
  axis equal;  %real ratios
  view(2);     %standard view in 2D
  hold off;
  axis off;
end