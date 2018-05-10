
function [coord,elem,surf,neumann,Q]=mesh_Q1(level,size_xy,size_z,size_hole)

% =========================================================================
%
%  This function creates hexahedral mesh for Q1 elements
%
%  input data:
%    level     - an integer defining a density of a uniform mesh
%    size_xy   - size of the body in directions x and y (integer)
%    size_z    - size of the body in z-direction (integer) 
%    size_hole - size of the hole in the body (integer)
%                size_hole < size_xy
%    body=(0,size_xy)  x(0,size_xy)  x(0,size_z)\setminus
%         (0,size_hole)x(0,size_hole)x(0,size_z)
%
%  output data:
%    coord   - coordinates of the nodes, size(coord)=(3,n_n) where n_n is a
%              number of nodes
%    elem    - array containing numbers of nodes defining each element,
%              size(elem)=(8,n_e), n_e = number of elements
%    surf    - array containing numbers of nodes defining each surace element,
%              size(surf)=(4,n_s), n_s = number of surface elements
%    neumann - array containing numbers of nodes defining each surface element,
%              size(neumann)=(4,n_e_s). The surface is the following side
%              of the body: (0,size_xy) x size_xy x (0,size_z), where the nonhomogeneous
%              Neumann boundary condition is considered.
%    Q       - logical array indicating the nodes where the homogeneous 
%              Dirichlet boundary condition is considered, size(Q)=(3,n_n)
%
% ======================================================================
%

%
% numbers of segments, nodes and elements
%

  N_x = size_xy*2^level;      % number of segments in x direction
  N_y = N_x;                  % number of segments in y direction
  N_z = size_z*2^level;       % number of segments in z direction
  N1_x=size_hole*2^level;     % number of segments in x direction ...    
  N2_x=N_x-N1_x;              % specifying the hole
  N1_y=size_hole*2^level;     % number of segments in y direction ...    
  N2_y=N_y-N1_y;              % specifying the hole
  % 
  n_node_xy = (N_x+1)*(N_y+1)-N1_x*N1_y;  % number of nodes in xy plane
  n_n = n_node_xy*(N_z+1);                % total number of nodes
 
%
% C - 3D auxilliary array that contains node numbers and that is important 
% for the mesh construction. Since the body is a union of two cuboids
% the array C also consists of two auxilliary 3D arrays, C1 and C2.
%
  C=zeros(N_x+1,N_y+1,N_z+1);
  C1=reshape(1:(N2_x+1)*N1_y*(N_z+1),N2_x+1,N1_y,N_z+1);
  C2=reshape(((N2_x+1)*N1_y*(N_z+1)+1):n_n,N_x+1,N2_y+1,N_z+1);
  C((N1_x+1):(N_x+1),1:N1_y,:)=C1;
  C(1:(N_x+1),(N1_y+1):(N_y+1),:)=C2;

%
% coordinates of nodes
%

  % coordinates in directions x, y and z
  coord_x=linspace(0,size_xy,N_x+1);
  coord_y=linspace(0,size_xy,N_y+1);
  coord_z=linspace(0,size_z,N_z+1);
  
  % long 1D arrays containing coordinates of all nodes in x,y,z directions
  c_x=[repmat(coord_x((N1_x+1):(N_x+1)),1,N1_y*(N_z+1)),...
           repmat(coord_x,1,(N2_y+1)*(N_z+1)) ];     
  c_y=[repmat(kron(coord_y(1:N1_y),ones(1,N2_x+1)),1,N_z+1),...
           repmat(kron(coord_y((N1_y+1):(N_y+1)),ones(1,N_x+1)),1,N_z+1)];     
  c_z=[kron(coord_z,ones(1,N1_y*(N2_x+1))),...
           kron(coord_z,ones(1,(N_x+1)*(N2_y+1)))]; 
       
  % the required array of coordinates, size(coord)=(3,n_n)
  coord=[c_x; c_y; c_z] ;
  
% 
% construction of the array elem
%

  % ordering of the nodes creating the unit cube:
  %  V1 -> [0 0 0], V2 -> [1 0 0], V3 -> [1 1 0], V4 -> [0 1 0]
  %  V5 -> [0 0 1], V6 -> [1 0 1], V7 -> [1 1 1], V8 -> [0 1 1]
  %  V1,...,V8 are logical 3D arrays which enable to select appropriate
  %  nodes from the array C.
  V1=false(N_x+1,N_y+1,N_z+1);
  V1((N1_x+1):N_x,1:N1_y,1:N_z)=1;
  V1(1:N_x,(N1_y+1):N_y,1:N_z)=1; 
  %
  V2=false(N_x+1,N_y+1,N_z+1);
  V2((N1_x+2):(N_x+1),1:N1_y,1:N_z)=1;
  V2(2:(N_x+1),(N1_y+1):N_y,1:N_z)=1;
  %
  V3=false(N_x+1,N_y+1,N_z+1);
  V3((N1_x+2):(N_x+1),2:(N1_y+1),1:N_z)=1;
  V3(2:(N_x+1),(N1_y+2):(N_y+1),1:N_z)=1;
  %
  V4=false(N_x+1,N_y+1,N_z+1);
  V4((N1_x+1):N_x,2:(N1_y+1),1:N_z)=1;
  V4(1:N_x,(N1_y+2):(N_y+1),1:N_z)=1;
  %
  V5=false(N_x+1,N_y+1,N_z+1);
  V5((N1_x+1):N_x,1:N1_y,2:(N_z+1))=1;
  V5(1:N_x,(N1_y+1):N_y,2:(N_z+1))=1;
  %
  V6=false(N_x+1,N_y+1,N_z+1);
  V6((N1_x+2):(N_x+1),1:N1_y,2:(N_z+1))=1;
  V6(2:(N_x+1),(N1_y+1):N_y,2:(N_z+1))=1;
  %
  V7=false(N_x+1,N_y+1,N_z+1);
  V7((N1_x+2):(N_x+1),2:(N1_y+1),2:(N_z+1))=1;
  V7(2:(N_x+1),(N1_y+2):(N_y+1),2:(N_z+1))=1;
  %
  V8=false(N_x+1,N_y+1,N_z+1);
  V8((N1_x+1):N_x,2:(N1_y+1),2:(N_z+1))=1;
  V8(1:N_x,(N1_y+2):(N_y+1),2:(N_z+1))=1;

  % the array elem, size(elem)=(8,n_e)
  elem=[C(V1)'; C(V2)'; C(V3)'; C(V4)';
        C(V5)'; C(V6)'; C(V7)'; C(V8)' ];
    
%
% Surface of the body - the array "surf"
%
  
  % For each face of the body, we define the restriction C_s of the array C 
  % and logical 2D arrays V1_s,...,V4_s which enable to select appropriate
  % nodes from the array C_s. We consider the following ordering of the 
  % nodes within the unit square:
  %   V1_s -> [0 0], V2_s -> [1 0], V3_s -> [1 1], V4_s -> [0 1]
  
  % Face 1: z=0 (the bottom of the body)
  C_s=zeros(N_x+1,N_y+1);
  C_s(:,:)=C(:,:,1);  
  V1_s=false(N_x+1,N_y+1);  V1_s((N1_x+1):N_x    ,1:N1_y    )=1; V1_s(1:N_x    ,(N1_y+1):N_y    )=1;
  V2_s=false(N_x+1,N_y+1);  V2_s((N1_x+2):(N_x+1),1:N1_y    )=1; V2_s(2:(N_x+1),(N1_y+1):N_y    )=1;
  V3_s=false(N_x+1,N_y+1);  V3_s((N1_x+2):(N_x+1),2:(N1_y+1))=1; V3_s(2:(N_x+1),(N1_y+2):(N_y+1))=1;
  V4_s=false(N_x+1,N_y+1);  V4_s((N1_x+1):N_x    ,2:(N1_y+1))=1; V4_s(1:N_x    ,(N1_y+2):(N_y+1))=1;
  surf1=[C_s(V1_s)'; C_s(V2_s)'; C_s(V3_s)'; C_s(V4_s)'];
  
  % Face 2: z=size_z (the top of the body)
  C_s=zeros(N_x+1,N_y+1);
  C_s(:,:)=C(:,:,end);  
  surf2=[C_s(V1_s)'; C_s(V2_s)'; C_s(V3_s)'; C_s(V4_s)'];
  
  % Face 3: y=0 (the front of the body)
  C_s=zeros(N_x+1,N_z+1);
  C_s(:,:)=C(:,1,:);  
  V1_s=false(N_x+1,N_z+1);  V1_s((N1_x+1):N_x    ,1:N_z    )=1;
  V2_s=false(N_x+1,N_z+1);  V2_s((N1_x+2):(N_x+1),1:N_z    )=1;
  V3_s=false(N_x+1,N_z+1);  V3_s((N1_x+2):(N_x+1),2:(N_z+1))=1;
  V4_s=false(N_x+1,N_z+1);  V4_s((N1_x+1):N_x    ,2:(N_z+1))=1;  
  surf3=[C_s(V1_s)'; C_s(V2_s)'; C_s(V3_s)'; C_s(V4_s)'];
  
  % Face 4: x=size_xy (the right hand side of the body)
  C_s=zeros(N_y+1,N_z+1);
  C_s(:,:)=C(end,:,:);  
  V1_s=false(N_y+1,N_z+1);  V1_s(1:N_y    ,1:N_z    )=1;
  V2_s=false(N_y+1,N_z+1);  V2_s(2:(N_y+1),1:N_z    )=1;
  V3_s=false(N_y+1,N_z+1);  V3_s(2:(N_y+1),2:(N_z+1))=1;
  V4_s=false(N_y+1,N_z+1);  V4_s(1:N_y    ,2:(N_z+1))=1;  
  surf4=[C_s(V1_s)'; C_s(V2_s)'; C_s(V3_s)'; C_s(V4_s)'];
  
  % Face 5: y=size_xy (the back of the body)
  C_s=zeros(N_x+1,N_z+1);
  C_s(:,:)=C(:,end,:);  
  V1_s=false(N_x+1,N_z+1);  V1_s(1:N_x    ,1:N_z)=1;
  V2_s=false(N_x+1,N_z+1);  V2_s(2:(N_x+1),1:N_z)=1;
  V3_s=false(N_x+1,N_z+1);  V3_s(2:(N_x+1),2:(N_z+1))=1;
  V4_s=false(N_x+1,N_z+1);  V4_s(1:N_x    ,2:(N_z+1))=1;  
  surf5=[C_s(V1_s)'; C_s(V2_s)'; C_s(V3_s)'; C_s(V4_s)'];
  
  % Face 6: x=0 (the left hand side of the body)
  C_s=zeros(N_y+1,N_z+1);
  C_s(:,:)=C(1,:,:);  
  V1_s=false(N_y+1,N_z+1);  V1_s((N1_y+1):N_y    ,1:N_z    )=1;
  V2_s=false(N_y+1,N_z+1);  V2_s((N1_y+2):(N_y+1),1:N_z    )=1;
  V3_s=false(N_y+1,N_z+1);  V3_s((N1_y+2):(N_y+1),2:(N_z+1))=1;
  V4_s=false(N_y+1,N_z+1);  V4_s((N1_y+1):N_y    ,2:(N_z+1))=1;  
  surf6=[C_s(V1_s)'; C_s(V2_s)'; C_s(V3_s)'; C_s(V4_s)'];
  
  % Face 7: y=size_hole (the hole face perpendicular to y-axis)
  C_s=zeros(N_x+1,N_z+1);
  C_s(:,:)=C(:,N1_y+1,:);  
  V1_s=false(N_x+1,N_z+1);  V1_s(1:N1_x    ,1:N_z    )=1;
  V2_s=false(N_x+1,N_z+1);  V2_s(2:(N1_x+1),1:N_z    )=1;
  V3_s=false(N_x+1,N_z+1);  V3_s(2:(N1_x+1),2:(N_z+1))=1;
  V4_s=false(N_x+1,N_z+1);  V4_s(1:N1_x    ,2:(N_z+1))=1;  
  surf7=[C_s(V1_s)'; C_s(V2_s)'; C_s(V3_s)'; C_s(V4_s)']; 
  
  % Face 8: x=size_hole (the hole face perpendicular to x-axis)
  C_s=zeros(N_y+1,N_z+1);
  C_s(:,:)=C(N1_x+1,:,:);  
  V1_s=false(N_y+1,N_z+1);  V1_s(1:N1_y    ,1:N_z    )=1;
  V2_s=false(N_y+1,N_z+1);  V2_s(2:(N1_y+1),1:N_z    )=1;
  V3_s=false(N_y+1,N_z+1);  V3_s(2:(N1_y+1),2:(N_z+1))=1;
  V4_s=false(N_y+1,N_z+1);  V4_s(1:N1_y    ,2:(N_z+1))=1;  
  surf8=[C_s(V1_s)'; C_s(V2_s)'; C_s(V3_s)'; C_s(V4_s)'];
  
  % the array "surf"
  surf = [surf1 surf2 surf3 surf4 surf5 surf6 surf7 surf8] ;    
 
%
% Boundary conditions
%
  
  % Nonhomogeneous Neumann boundary conditions on Face 5
  neumann=surf5;     

  % logical array indicating the nodes with the Dirichlet boundary cond.
  Q = coord>0 ;
  Q(3,(coord(3,:)==size_z)) = 0;     
  
end