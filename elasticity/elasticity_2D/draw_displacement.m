function draw_displacement(coord,elem,U,U_disp,elem_type,...
                                    size_xy,size_hole)

% =========================================================================
%
%  This function depicts prescribed displacements, deformed and undeformed
%  shape of the body
%
%  input data:
%    coord     - coordinates of the nodes, size(coord)=(2,n_n) where n_n is a
%                number of nodes
%    elem      - array containing numbers of nodes defining each element,
%                size(elem)=(n_p,n_e), n_e = number of elements
%    U         - nodal displacements, size(U)=(2,n_n)
%    U_disp    - prescribed displacements (e.g. total displacement or 
%                displacements in x directions, etc.), size(U_disp)=(1,n_n)
%    elem_type - the type of finite elements; available choices:
%                'P1', 'P2', 'Q1', 'Q2'
%    size_xy   - size of the body in x and y direction (integer)
%    size_hole - size of the hole in the body (integer)
%                size_hole < size_xy
%    body=(0,size_xy)  x(0,size_xy)  \setminus
%         (0,size_hole)x(0,size_hole)
%
% ======================================================================
%

  figure;
  hold on;
  coord_aux = [coord;zeros(1,size(coord,2))] + [U;zeros(1,size(U,2))];
  % displacements and deformed shape
  if strcmp(elem_type,'P1')||strcmp(elem_type,'P2')
    s = patch('Faces',elem(1:3,:)','Vertices',coord_aux',...
        'FaceVertexCData',U_disp','FaceColor','interp','EdgeColor','none'); 
  else
    s = patch('Faces',elem(1:4,:)','Vertices',coord_aux',...
        'FaceVertexCData',U_disp','FaceColor','interp','EdgeColor','none'); 
  end
  if ~isOctave  
    alpha(s,.5);
  end
  colorbar;
  
  % undeformed shape of the body
  plot([size_hole,size_xy],[0,0])
  plot([0,size_xy],[size_xy,size_xy])
  plot([0,size_hole],[size_hole,size_hole])
  plot([0,0],[size_hole,size_xy])
  plot([size_hole,size_hole],[0,size_hole])
  plot([size_xy,size_xy],[0,size_xy]) 
  
  %
  box on
  view(2);      % standard view ve 2D
  axis equal;   % real ratios
  hold off;
  axis off;
end