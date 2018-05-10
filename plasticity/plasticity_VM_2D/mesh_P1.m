
function [coord,elem,surf,neumann,Q]=mesh_P1(level,size_xy,size_hole)

% =========================================================================
%
%  This function creates triangular mesh for P1 elements
%
%  input data:
%    level     - an integer defining a density of a uniform mesh
%    size_xy   - size of the body in directions x and y (integer)
%    size_hole - size of the hole in the body (integer)
%                size_hole < size_xy
%    body=(0,size_xy)  x(0,size_xy)  \setminus
%         (0,size_hole)x(0,size_hole)
%
%  output data:
%    coord   - coordinates of the nodes, size(coord)=(2,n_n) where n_n is a
%              number of nodes
%    elem    - array containing numbers of nodes defining each element, 
%              size(elem)=(3,n_e), n_e = number of elements
%    surf    - array containing numbers of nodes defining each surface element, 
%              size(surf)=(2,n_s), n_s = number of surface elements
%    neumann - array containing numbers of nodes defining each surface element, 
%              size(neuman)=(2,n_e_s). The surface is the following side 
%              of the body: (0,size_xy) x size_xy, where the nonhomogeneous
%              Neumann boundary condition is considered.
%    Q       - logical array indicating the nodes where the Dirichlet
%              boundary condition is considered, size(Q)=(2,n_n)
%
% ======================================================================
%

%
% numbers of segments, nodes and elements
%

  N_x = size_xy*2^level;      % number of segments in x direction
  N_y = N_x;                  % number of segments in y direction
  N1_x=size_hole*2^level;     % number of segments in x direction ...    
  N2_x=N_x-N1_x;              % specifying the hole
  N1_y=size_hole*2^level;     % number of segments in y direction ...    
  N2_y=N_y-N1_y;              % specifying the hole
  % 
  n_n = (N_x+1)*(N_y+1)-N1_x*N1_y;        % total number of nodes
  n_cell_xy = N_x*N_y-N1_x*N1_y;          % number of cells in xy plane
  n_e = n_cell_xy*2    ;                  % total number of elements

%
% C - 2D auxilliary array that contains node numbers and that is important 
% for the mesh construction. Since the body is a union of two rectangles
% the array C also consists of two auxilliary 2D arrays, C1 and C2.
%
  C=zeros(N_x+1,N_y+1);
  C1=reshape(1:(N2_x+1)*N1_y,N2_x+1,N1_y);
  C2=reshape(((N2_x+1)*N1_y+1):n_n,N_x+1,N2_y+1);
  C((N1_x+1):(N_x+1),1:N1_y)=C1;
  C(1:(N_x+1),(N1_y+1):(N_y+1))=C2;
  
%
% coordinates of nodes
%
  % coordinates in directions x and y
  coord_x=linspace(0,size_xy,N_x+1);
  coord_y=linspace(0,size_xy,N_y+1);
  
  % long 1D arrays containing coordinates of all nodes in x,y directions
  c_x=[repmat(coord_x((N1_x+1):(N_x+1)),1,N1_y),...
           repmat(coord_x,1,(N2_y+1))  ];     
  c_y=[repmat(kron(coord_y(1:N1_y),ones(1,N2_x+1)),1),...
           repmat(kron(coord_y((N1_y+1):(N_y+1)),ones(1,N_x+1)),1)];  
         
  % the required array of coordinates, size(coord)=(2,n_n)   
  coord=[c_x; c_y] ;  

% 
% construction of the array elem
%
  % ordering of the nodes creating the unit cube:
  %  V1 -> [0 0], V2 -> [1 0], V3 -> [1 1], V4 -> [0 1]
  %  V1,...,V4 are logical 2D arrays which enable to select appropriate
  %  nodes from the array C.

  V1=false(N_x+1,N_y+1);
  V1((N1_x+1):N_x,1:N1_y)=1;
  V1(1:N_x,(N1_y+1):N_y) =1; 
  %
  V2=false(N_x+1,N_y+1);
  V2((N1_x+2):(N_x+1),1:N1_y)=1;
  V2(2:(N_x+1),(N1_y+1):N_y) =1;
  %
  V3=false(N_x+1,N_y+1);
  V3((N1_x+2):(N_x+1),2:(N1_y+1))=1;
  V3(2:(N_x+1),(N1_y+2):(N_y+1)) =1;
  %
  V4=false(N_x+1,N_y+1);
  V4((N1_x+1):N_x,2:(N1_y+1))=1;
  V4(1:N_x,(N1_y+2):(N_y+1)) =1;
 

  % used division of a rectangle into 2 triangles:   
  %   V1 V2 V4
  %   V2 V3 V4
  % size(aux_elem)=(2*3,n_e/2)
  aux_elem=[C(V1)'; C(V2)'; C(V4)';
            C(V2)'; C(V3)'; C(V4)' ];
        
  % the array elem, size(elem)=(3,n_e)          
  elem=reshape(aux_elem,3,n_e);     
 
%
% Surface of the body - the array "surf"
%
  
  % For each face of the body, we define the restriction C_s of the array C 
  % and logical 2D arrays V1_s, V2_s which enable to select appropriate
  % nodes from the array C_s. We consider the following ordering of the 
  % nodes within the unit line:
  %   V1_s -> [0 0], V2_s -> [1 0]      
  
  % Face 1: y=0 (the front of the body)
  C_s=zeros(N_x+1,1);
  C_s(:)=C(:,1);  
  V1_s=false(N_x+1,1);  V1_s((N1_x+1):N_x    ,1)=1;
  V2_s=false(N_x+1,1);  V2_s((N1_x+2):(N_x+1),1)=1;  
  aux_surf=[C_s(V1_s)'; C_s(V2_s)' ];
  surf1=reshape(aux_surf,2,N2_x);       
  
  % Face 2: x=size_xy (the right hand side of the body)
  C_s=zeros(N_y+1,1);
  C_s(:)=C(end,:);  
  V1_s=false(N_y+1,1);  V1_s(1:N_y    ,1)=1;
  V2_s=false(N_y+1,1);  V2_s(2:(N_y+1),1)=1;
  aux_surf=[C_s(V1_s)'; C_s(V2_s)' ];
  surf2=reshape(aux_surf,2,N_y);       
  
  % Face 3: y=size_xy (the back of the body)
  C_s=zeros(N_x+1,1);
  C_s(:)=C(:,end);  
  V1_s=false(N_x+1,1);  V1_s(1:N_x    ,1)=1;
  V2_s=false(N_x+1,1);  V2_s(2:(N_x+1),1)=1;  
  aux_surf=[C_s(V1_s)'; C_s(V2_s)' ];
  surf3=reshape(aux_surf,2,N_x);       
  
  % Face 4: x=0 (the left hand side of the body)
  C_s=zeros(N_y+1,1);
  C_s(:)=C(1,:);  
  V1_s=false(N_y+1,1);  V1_s((N1_y+1):N_y    ,1)=1;
  V2_s=false(N_y+1,1);  V2_s((N1_y+2):(N_y+1),1)=1;  
  aux_surf=[C_s(V1_s)'; C_s(V2_s)' ];
  surf4=reshape(aux_surf,2,N2_y);     
  
  % Face 5: y=size_hole (the hole face perpendicular to y-axis)
  C_s=zeros(N_x+1,1);
  C_s(:)=C(:,N1_y+1);  
  V1_s=false(N_x+1,1);  V1_s(1:N1_x    ,1)=1;
  V2_s=false(N_x+1,1);  V2_s(2:(N1_x+1),1)=1; 
  aux_surf=[C_s(V1_s)'; C_s(V2_s)' ];
  surf5=reshape(aux_surf,2,N1_x);       
  
  % Face 6: x=size_hole (the hole face perpendicular to x-axis)
  C_s=zeros(N_y+1,1);
  C_s(:)=C(N1_x+1,:);  
  V1_s=false(N_y+1,1);  V1_s(1:N1_y    ,1)=1;
  V2_s=false(N_y+1,1);  V2_s(2:(N1_y+1),1)=1;
  aux_surf=[C_s(V1_s)'; C_s(V2_s)' ];
  surf6=reshape(aux_surf,2,N1_y);  
  
  % the array "surf"
  surf = [surf1 surf2 surf3 surf4 surf5 surf6] ;
  
%
% Boundary conditions
%
  
  % nonhomogeneous Neumann boundary conditions on Face 3
  neumann=surf3;          

  % logical array indicating the nodes with the Dirichlet boundary cond.
  Q = coord>0 ;
  
end
