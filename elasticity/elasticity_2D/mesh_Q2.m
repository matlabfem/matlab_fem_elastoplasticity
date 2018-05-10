
function [coord,elem,surf,neumann,dirichlet,Q]=...
                    mesh_Q2(level,size_xy,size_hole)

% =========================================================================
%
%  This function creates quadrilateral mesh for Q2 elements consisting of 
%  4 vertices and 4 midpoints.
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
%    coord     - coordinates of the nodes, size(coord)=(2,n_n) where n_n
%                is a number of nodes including midpoints
%    elem      - array containing numbers of nodes defining each element,
%                size(elem)=(8,n_e), n_e = number of elements
%    surf      - array containing numbers of nodes defining each surface
%                element, size(surf)=(3,n_s), n_s = number of surface elements
%    neumann   - array containing numbers of nodes defining each
%                surface element, size(neuman)=(3,n_e_s). The surface 
%                is the following side of the body: (0,size_xy) x size_xy, 
%                where the nonhomogeneous Neumann boundary condition 
%                is considered.
%    dirichlet - array indicating the nodes where the
%                nonhomogeneous Dirichlet boundary condition is considered,
%                size(dirichlet)=(2,n_n)
%    Q         - logical array indicating the nodes where the Dirichlet
%                boundary condition is considered, size(Q)=(2,n_n)
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
% C - 2D auxilliary array that contains node numbers and that is important 
% for the mesh construction. Since the body is a union of two 
% the array C also consists of two auxilliary 2D arrays, C1 and C2. To
% specify selected midpoints, we define 2D logical arrays Q1 and Q2.
%
  C=zeros(2*N_x+1,2*N_y+1);
  %
  C1=zeros(2*N2_x+1,2*N1_y);
  Q1= true(2*N2_x+1,2*N1_y);
  Q1(2:2:end,2:2:end)=0;
  C1(Q1)=1:length(C1(Q1));
  %   
  C2=zeros(2*N_x+1,2*N2_y+1);
  Q2= true(2*N_x+1,2*N2_y+1);
  Q2(2:2:end,2:2:end)=0;
  C2(Q2)=length(C1(Q1))+(1:length(C2(Q2)));
  %
  C((2*N1_x+1):(2*N_x+1),1:2*N1_y)=C1;
  C(1:(2*N_x+1),(2*N1_y+1):(2*N_y+1))=C2;
  
%
% coordinates of nodes
%
  % coordinates in directions x, y and z
  coord_x=linspace(0,size_xy,2*N_x+1);
  coord_y=linspace(0,size_xy,2*N_y+1);
  
  % 2D arrays containing coordinates in x and y directions
  c1_x=repmat(coord_x((2*N1_x+1):(2*N_x+1))',[1,2*N1_y]);
  c2_x=repmat(coord_x',[1,(2*N2_y+1)]);
  %
  c1_y=repmat(coord_y(1:2*N1_y),[2*N2_x+1,1]);
  c2_y=repmat(coord_y((2*N1_y+1):(2*N_y+1)),[2*N_x+1,1]);  
  
  % the required array of coordinates, size(coord)=(2,n_n)
  coord=[c1_x(Q1)' c2_x(Q2)'; 
         c1_y(Q1)' c2_y(Q2)' ];
  
% 
% construction of the array elem, size(elem)=(8,n_e)
%
  % ordering of the nodes creating the unit cube:
  %  V1 -> [0 0], V2 -> [1 0], V3 -> [1 1], V4 -> [0 1]
  %  V1,...,V4 are logical 2D arrays which enable to select appropriate
  %  nodes from the array C.

  V1=false(2*N_x+1,2*N_y+1);
  V1((2*N1_x+1):2:(2*N_x-1),         1:2:(2*N1_y-1))=1;
  V1(         1:2:(2*N_x-1),(2*N1_y+1):2:(2*N_y-1) )=1;
  %
  V2=false(2*N_x+1,2*N_y+1);
  V2((2*N1_x+3):2:(2*N_x+1),         1:2:(2*N1_y-1))=1;
  V2(         3:2:(2*N_x+1),(2*N1_y+1):2:(2*N_y-1) )=1;
  %
  V3=false(2*N_x+1,2*N_y+1);
  V3((2*N1_x+3):2:(2*N_x+1),         3:2:(2*N1_y+1))=1;
  V3(         3:2:(2*N_x+1),(2*N1_y+3):2:(2*N_y+1) )=1;
  %
  V4=false(2*N_x+1,2*N_y+1);
  V4((2*N1_x+1):2:(2*N_x-1),         3:2:(2*N1_y+1))=1;
  V4(         1:2:(2*N_x-1),(2*N1_y+3):2:(2*N_y+1) )=1;
  
  % logical arrays for midpoints, e.g. V12 represents the midpoints between
  % V1 and V2
  V12=false(2*N_x+1,2*N_y+1);
  V12((2*N1_x+2):2:(2*N_x),         1:2:(2*N1_y-1))=1;
  V12(         2:2:(2*N_x),(2*N1_y+1):2:(2*N_y-1) )=1;
  %
  V14=false(2*N_x+1,2*N_y+1);
  V14((2*N1_x+1):2:(2*N_x-1),         2:2:(2*N1_y))=1;
  V14(         1:2:(2*N_x-1),(2*N1_y+2):2:(2*N_y) )=1;
  %
  V23=false(2*N_x+1,2*N_y+1);
  V23((2*N1_x+3):2:(2*N_x+1),         2:2:(2*N1_y))=1;
  V23(         3:2:(2*N_x+1),(2*N1_y+2):2:(2*N_y) )=1;
  %
  V34=false(2*N_x+1,2*N_y+1);
  V34((2*N1_x+2):2:(2*N_x),         3:2:(2*N1_y+1))=1;
  V34(         2:2:(2*N_x),(2*N1_y+3):2:(2*N_y+1) )=1;

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
  
  % Face 1: y=0 (the front of the body)
  C_s=zeros(2*N_x+1,1);
  C_s(:)=C(:,1);
   V1_s=false(2*N_x+1,1);   V1_s((2*N1_x+1):2:(2*N_x-1),1)=1;
   V2_s=false(2*N_x+1,1);   V2_s((2*N1_x+3):2:(2*N_x+1),1)=1;
  V12_s=false(2*N_x+1,1);  V12_s((2*N1_x+2):2:(2*N_x  ),1)=1;
  surf1=[ C_s(V1_s)';  C_s(V2_s)'; C_s(V12_s)' ];
  
  % Face 2: x=size_xy (the right hand side of the body)
  C_s=zeros(2*N_y+1,1);
  C_s(:)=C(end,:);  
   V1_s=false(2*N_y+1,1);   V1_s(1:2:(2*N_y-1),1)=1;
   V2_s=false(2*N_y+1,1);   V2_s(3:2:(2*N_y+1),1)=1;
  V12_s=false(2*N_y+1,1);  V12_s(2:2:(2*N_y  ),1)=1;
  surf2=[ C_s(V1_s)';  C_s(V2_s)';  C_s(V12_s)' ];
  
  % Face 3: y=size_xy (the back of the body)
  C_s=zeros(2*N_x+1,1);
  C_s(:)=C(:,end);
   V1_s=false(2*N_x+1,1);   V1_s(1:2:(2*N_x-1),1)=1;
   V2_s=false(2*N_x+1,1);   V2_s(3:2:(2*N_x+1),1)=1;
  V12_s=false(2*N_x+1,1);  V12_s(2:2:(2*N_x  ),1)=1;
  surf3=[ C_s(V1_s)';  C_s(V2_s)';  C_s(V12_s)' ];
  
  % Face 4: x=0 (the left hand side of the body)
  C_s=zeros(2*N_y+1,1);
  C_s(:)=C(1,:);  
   V1_s=false(2*N_y+1,1);   V1_s((2*N1_y+1):2:(2*N_y-1),1)=1;
   V2_s=false(2*N_y+1,1);   V2_s((2*N1_y+3):2:(2*N_y+1),1)=1;
  V12_s=false(2*N_y+1,1);  V12_s((2*N1_y+2):2:(2*N_y  ),1)=1;
  surf4=[ C_s(V1_s)';  C_s(V2_s)';  C_s(V12_s)' ];
  
  % Face 5: y=size_hole (the hole face perpendicular to y-axis)
  C_s=zeros(2*N_x+1,1);
  C_s(:)=C(:,2*N1_y+1);
   V1_s=false(2*N_x+1,1);   V1_s(1:2:(2*N1_x-1),1)=1;
   V2_s=false(2*N_x+1,1);   V2_s(3:2:(2*N1_x+1),1)=1;
  V12_s=false(2*N_x+1,1);  V12_s(2:2:(2*N1_x  ),1)=1;
  surf5=[ C_s(V1_s)';  C_s(V2_s)';  C_s(V12_s)' ];        
  
  % Face 6: x=size_hole (the hole face perpendicular to x-axis)
  C_s=zeros(2*N_y+1,1);
  C_s(:)=C(2*N1_x+1,:);  
   V1_s=false(2*N_y+1,1);   V1_s(1:2:(2*N1_y-1),1)=1;
   V2_s=false(2*N_y+1,1);   V2_s(3:2:(2*N1_y+1),1)=1;
  V12_s=false(2*N_y+1,1);  V12_s(2:2:(2*N1_y  ),1)=1;
  surf6=[ C_s(V1_s)';  C_s(V2_s)';  C_s(V12_s)' ];  
  
  % the array "surf"
  surf = [surf1 surf2 surf3 surf4 surf5 surf6] ;    
  
%
% Boundary conditions
%
  
  % Nonhomogeneous Neumann boundary conditions on Face 3
  neumann=surf3;       
   
  % array indicating the nodes with non-homogen. Dirichlet boundary cond.
  dirichlet = zeros(size(coord));
  dirichlet(1,(coord(2,:) == 0)) = 1;  

  % logical array indicating the nodes with the Dirichlet boundary cond.
  Q = coord>0 ;
  Q(1,(coord(2,:) == 0)) = 0;
  
end
