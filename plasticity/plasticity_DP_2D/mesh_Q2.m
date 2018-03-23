
function [coord,elem,surf,dirichlet,Q]=mesh_Q2(level,size_xy)

% =========================================================================
%
%  This function creates quadrilateral mesh for Q2 elements consisting of 
%  4 vertices and 4 midpoints.
%
%  input data:
%    level - an integer defining a density of a uniform mesh
%    size_xy - size of the body in directions x and y (integer)
%         body=(0,size_xy)x(0,size_xy)
%
%  output data:
%    coord - coordinates of the nodes, size(coord)=(2,n_n) where n_n is a
%            number of nodes including midpoints
%    elem - 8 x n_e array containing numbers of nodes defining each
%           element, n_e = number of elements
%    surf - 3 x n_s array containing numbers of nodes defining each
%           surface element, n_s = number of surface elements
%    dirichlet - 2 x n_n array indicating the nodes where the
%            nonhomogeneous Dirichlet boundary condition is considered
%    Q - logical 2 x n_n array indicating the nodes where the Dirichlet
%        boundary condition is considered
%
% ======================================================================
%

%
% numbers of segments, nodes and elements
%

  N_x = size_xy*2^level;      % number of segments in x direction
  N_y = N_x;                  % number of segments in y direction

%
% C - 2D auxilliary array that contains node numbers and that is important 
% for the mesh construction. Since the body is a union of two 
% the array C also consists of two auxilliary 2D arrays, C1 and C2. To
% specify selected midpoints, we define 2D logical arrays Q1 and Q2.
%
  C=zeros(2*N_x+1,2*N_y+1);
  Q_mid= true(2*N_x+1,2*N_y+1);
  Q_mid(2:2:end,2:2:end)=0;
  C(Q_mid)=1:length(C(Q_mid));
  
%
% coordinates of nodes
%
  % coordinates in directions x, y and z
  coord_x=linspace(0,size_xy,2*N_x+1);
  coord_y=linspace(0,size_xy,2*N_y+1);
  
  % 2D arrays containing coordinates in x and y directions
  c_x=repmat(coord_x',[1,2*N_y+1]);
  %
  c_y=repmat(coord_y,[2*N_x+1,1]);
  
  % the required 3 x n_n array of coordinates
  coord=[c_x(Q_mid)'; c_y(Q_mid)'];
  
% 
% construction of the 8 x n_e array elem
%
  % ordering of the nodes creating the unit cube:
  %  V1 -> [0 0], V2 -> [1 0], V3 -> [1 1], V4 -> [0 1]
  %  V1,...,V4 are logical 2D arrays which enable to select appropriate
  %  nodes from the array C.

  V1=false(2*N_x+1,2*N_y+1);
  V1(1:2:(2*N_x-1),1:2:(2*N_y-1))=1;
  %
  V2=false(2*N_x+1,2*N_y+1);
  V2(3:2:(2*N_x+1),1:2:(2*N_y-1))=1;
  %
  V3=false(2*N_x+1,2*N_y+1);
  V3(3:2:(2*N_x+1),3:2:(2*N_y+1))=1;
  %
  V4=false(2*N_x+1,2*N_y+1);
  V4(1:2:(2*N_x-1),3:2:(2*N_y+1))=1;
  
  % logical arrays for midpoints, e.g. V12 represents the midpoints between
  % V1 and V2
  V12=false(2*N_x+1,2*N_y+1);
  V12(2:2:(2*N_x),1:2:(2*N_y-1))=1;
  %
  V14=false(2*N_x+1,2*N_y+1);
  V14(1:2:(2*N_x-1),2:2:(2*N_y))=1;
  %
  V23=false(2*N_x+1,2*N_y+1);
  V23(3:2:(2*N_x+1),2:2:(2*N_y))=1;
  %
  V34=false(2*N_x+1,2*N_y+1);
  V34(2:2:(2*N_x),3:2:(2*N_y+1))=1;
  
  % used ordering of nodes within a cube:   
  % V1 V2 V3 V4 V12 V23 V34 V14
  
  elem=[ C(V1)';  C(V2)';  C(V3)';  C(V4)';  
        C(V12)'; C(V23)'; C(V34)'; C(V14)' ];
    
%
% Surface of the body - the array "surf"
%
  
  % For each face of the body, we define the restriction C_s of the array C 
  % and logical 2D arrays V1_s,V2_s,V12_s which enable to select appropriate
  % nodes from the array C_s. We consider the following ordering of the 
  % nodes within the unit square:
  %      V1_s -> [0 0], V2_s -> [1 0], V12_s -> [1/2 0]
  
 % Edge 1: y=0 (the bottom of the body)
  C_s=zeros(2*N_x+1,1);
  C_s(:)=C(:,1);
   V1_s=false(2*N_x+1,1);   V1_s(1:2:(2*N_x-1),1)=1;
   V2_s=false(2*N_x+1,1);   V2_s(3:2:(2*N_x+1),1)=1;
  V12_s=false(2*N_x+1,1);  V12_s(2:2:(2*N_x  ),1)=1;
  aux_surf=[C_s(V1_s)'; C_s(V2_s)'; C_s(V12_s)' ];
  surf1=reshape(aux_surf,3,N_x);     
  
  % Edge 2: x=size_xy (the right hand side of the body)
  C_s=zeros(2*N_y+1,1);
  C_s(:)=C(end,:);  
   V1_s=false(2*N_y+1,1);   V1_s(1:2:(2*N_y-1),1)=1;
   V2_s=false(2*N_y+1,1);   V2_s(3:2:(2*N_y+1),1)=1;
  V12_s=false(2*N_y+1,1);  V12_s(2:2:(2*N_y  ),1)=1;
  aux_surf=[C_s(V1_s)'; C_s(V2_s)'; C_s(V12_s)' ];
  surf2=reshape(aux_surf,3,N_y);       
  
  % Edge 3: y=size_xy (the top of the body)
  C_s=zeros(2*N_x+1,1);
  C_s(:)=C(:,end);
   V1_s=false(2*N_x+1,1);   V1_s(1:2:(2*N_x-1),1)=1;
   V2_s=false(2*N_x+1,1);   V2_s(3:2:(2*N_x+1),1)=1;
  V12_s=false(2*N_x+1,1);  V12_s(2:2:(2*N_x  ),1)=1;
  aux_surf=[C_s(V1_s)'; C_s(V2_s)'; C_s(V12_s)' ];
  surf3=reshape(aux_surf,3,N_x);     
  
  % Edge 4: x=0 (the left hand side of the body)
  C_s=zeros(2*N_y+1,1);
  C_s(:)=C(1,:);  
   V1_s=false(2*N_y+1,1);   V1_s(1:2:(2*N_y-1),1)=1;
   V2_s=false(2*N_y+1,1);   V2_s(3:2:(2*N_y+1),1)=1;
  V12_s=false(2*N_y+1,1);  V12_s(2:2:(2*N_y  ),1)=1;
  aux_surf=[C_s(V1_s)'; C_s(V2_s)'; C_s(V12_s)' ];
  surf4=reshape(aux_surf,3,N_y);       
  
  % the array "surf"
  surf = [surf1 surf2 surf3 surf4] ;
  
%
% Boundary conditions
%   

  % array indicating the nodes with non-homogen. Dirichlet boundary cond.
  dirichlet = zeros(size(coord));
  dirichlet(2,(coord(2,:)==size_xy)&(coord(1,:)<=1.0001)) = 1;     
  
  % logical array indicating the nodes with the Dirichlet boundary cond.
  Q = coord>0 ;
  Q(2,(coord(2,:)==size_xy)&(coord(1,:)<=1.0001)) = 0;      
  Q(1,(coord(1,:)==size_xy)) = 0;   
 
end
