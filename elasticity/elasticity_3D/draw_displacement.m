function draw_displacement(coord,surf,U,U_disp,elem_type,...
                                    size_xy,size_z,size_hole)

% =========================================================================
%
%  This function depicts prescribed displacements, deformed and undeformed
%  shape of the body
%
%  input data:
%    coord     - coordinates of the nodes, size(coord)=(3,n_n) where n_n is a
%                number of nodes
%    surf      - array containing numbers of nodes defining eachsurface element,
%                size(surf)=(n_p,n_s), n_s = number of surface elements
%    U         - nodal displacements, size(U)=(3,n_n)
%    U_disp    - prescribed displacements (e.g. total displacement or 
%                displacements in x directions, etc.), size(U_disp)=(1,n_n)
%    elem_type - the type of finite elements; available choices:
%                'P1', 'P2', 'Q1', 'Q2'
%    size_xy   - size of the body in x and y direction (integer)
%    size_z    - size of the body in z-direction (integer) 
%    size_hole - size of the hole in the body (integer)
%                size_hole < size_xy
%    body=(0,size_xy)  x(0,size_xy)  x(0,size_z)\setminus
%         (0,size_hole)x(0,size_hole)x(0,size_z)
%
% ======================================================================
%

  figure;
  hold on;
  
  % displacements and deformed shape
  if strcmp(elem_type,'P1')||strcmp(elem_type,'P2')
     s = patch('Faces',surf(1:3,:)','Vertices',coord'+U',...
        'FaceVertexCData',U_disp','FaceColor','interp','EdgeColor','none'); 
  else
     s = patch('Faces',surf(1:4,:)','Vertices',coord'+U',...
        'FaceVertexCData',U_disp','FaceColor','interp','EdgeColor','none'); 
  end
  if ~isOctave  
    alpha(s,.5);
  end
  colorbar;
  
  % undeformed shape of the body
  plot3([size_hole,size_xy],[0,0],[0,0])
  plot3([size_hole,size_xy],[0,0],[size_z,size_z])
  plot3([0,size_xy],[size_xy,size_xy],[0,0])
  plot3([0,size_xy],[size_xy,size_xy],[size_z,size_z])
  plot3([0,size_hole],[size_hole,size_hole],[0,0])
  plot3([0,size_hole],[size_hole,size_hole],[size_z,size_z])
  plot3([0,0],[size_hole,size_xy],[0,0])
  plot3([0,0],[size_hole,size_xy],[size_z,size_z])
  plot3([size_hole,size_hole],[0,size_hole],[0,0])
  plot3([size_hole,size_hole],[0,size_hole],[size_z,size_z])
  plot3([size_xy,size_xy],[0,size_xy],[0,0])
  plot3([size_xy,size_xy],[0,size_xy],[size_z,size_z])
  plot3([0,0],[size_hole,size_hole],[0,size_z])
  plot3([0,0],[size_xy,size_xy],[0,size_z])
  plot3([size_hole,size_hole],[0,0],[0,size_z])
  plot3([size_hole,size_hole],[size_hole,size_hole],[0,size_z])  
  plot3([size_xy,size_xy],[0,0],[0,size_z])
  plot3([size_xy,size_xy],[size_xy,size_xy],[0,size_z])  
  
  %
  box on
  view(3);      % standard view ve 3D
  axis equal;   % real ratios
  hold off;
  axis off;
end