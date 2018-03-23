
function [coord,elem,surf,neumann,Q]=mesh_P2(level,size_xy,size_z,size_hole)

% =========================================================================
%
%  This function creates tetrahedral mesh for P2 elements
%
%  input data:
%    level - an integer defining a density of a uniform mesh
%    size_xy - size of the body in directions x and y (integer)
%    size_z  - size of the body in z-direction (integer) 
%    size_hole - size of the hole in the body (integer)
%                size_hole < size_xy
%         body=(0,size_xy)x(0,size_xy)x(0,size_z)\setminus
%              (0,size_hole)x(0,size_hole)x(0,size_z)
%
%  output data:
%    coord - coordinates of the nodes, size(coord)=(3,n_n) where n_n is a
%            number of nodes including midpoints
%    elem - 10 x n_e array containing numbers of nodes defining each
%           element, n_e = number of elements
%    surf - 6 x n_s array containing numbers of nodes defining each
%           surface element, n_s = number of surface elements
%    neumann - 6 x n_e_s array containing numbers of nodes defining each
%           surface element. The surface is the following side of the body:
%           (0,size_xy) x size_xy x (0,size_z), where the nonhomogeneous
%           Neumann boundary condition is considered.
%    Q - logical 3 x n_n array indicating the nodes where the homogeneous 
%        Dirichlet boundary condition is considered
%
% ======================================================================
%

%
% numbers of segments, nodes and elements
%

  N_x = size_xy*2^level;      % number of segments in x direction
  N_y = N_x;                  % number of segments in y direction
  N_z = size_z*2^level;       % number of segments in z-direction
  N1_x=size_hole*2^level;     % number of segments in x direction ...    
  N2_x=N_x-N1_x;              % specifying the hole
  N1_y=size_hole*2^level;     % number of segments in y direction ...    
  N2_y=N_y-N1_y;              % specifying the hole
  % 
  n_node_xy=(2*N_x+1)*(2*N_y+1)-4*N1_x*N1_y; % number of nodes in xy plane
  n_n = n_node_xy*(2*N_z+1);                 % total number of nodes
  n_cell_xy = N_x*N_y-N1_x*N1_y;             % number of cells in xy plane
  n_e = n_cell_xy*N_z*6;                     % total number of elements
  
%
% C - 3D auxilliary array that contains node numbers and that is important 
% for the mesh construction. Since the body is a union of two cuboids
% the array C also consists of two auxilliary 3D arrays, C1 and C2.
%
  C=zeros(2*N_x+1,2*N_y+1,2*N_z+1);
  C1=reshape(1:(2*N2_x+1)*2*N1_y*(2*N_z+1),2*N2_x+1,2*N1_y,2*N_z+1);
  C2=reshape(((2*N2_x+1)*2*N1_y*(2*N_z+1)+1):n_n,2*N_x+1,2*N2_y+1,2*N_z+1);
  C((2*N1_x+1):(2*N_x+1),1:2*N1_y,:)=C1;
  C(1:(2*N_x+1),(2*N1_y+1):(2*N_y+1),:)=C2;
  
%
% coordinates of nodes
%
  % coordinates in directions x, y and z
  coord_x=linspace(0,size_xy,2*N_x+1);
  coord_y=linspace(0,size_xy,2*N_y+1);
  coord_z=linspace(0,size_z,2*N_z+1);
  
  % long 1D arrays containing coordinates of all nodes in x,y,z directions
  c_x=[repmat(coord_x((2*N1_x+1):(2*N_x+1)),1,2*N1_y*(2*N_z+1)),...
       repmat(coord_x,1,(2*N2_y+1)*(2*N_z+1)) ];     
  c_y=[repmat(kron(coord_y(1:2*N1_y),ones(1,2*N2_x+1)),1,2*N_z+1),...
       repmat(kron(coord_y((2*N1_y+1):(2*N_y+1)),ones(1,2*N_x+1)),1,2*N_z+1)];     
  c_z=[kron(coord_z,ones(1,2*N1_y*(2*N2_x+1))),...
       kron(coord_z,ones(1,(2*N_x+1)*(2*N2_y+1)))];    
       
  % the required 3 x n_n array of coordinates          
  coord=[c_x; c_y; c_z] ;
  
% 
% construction of the array elem
%
  % ordering of the nodes creating the unit cube:
  %  V1 -> [0 0 0], V2 -> [1 0 0], V3 -> [1 1 0], V4 -> [0 1 0]
  %  V5 -> [0 0 1], V6 -> [1 0 1], V7 -> [1 1 1], V8 -> [0 1 1]
  %  V1,...,V8 are logical 3D arrays which enable to select appropriate
  %  nodes from the array C.

  V1=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V1((2*N1_x+1):2:(2*N_x-1),         1:2:(2*N1_y-1),1:2:(2*N_z-1))=1;
  V1(         1:2:(2*N_x-1),(2*N1_y+1):2:(2*N_y-1) ,1:2:(2*N_z-1))=1;
  %
  V2=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V2((2*N1_x+3):2:(2*N_x+1),         1:2:(2*N1_y-1),1:2:(2*N_z-1))=1;
  V2(         3:2:(2*N_x+1),(2*N1_y+1):2:(2*N_y-1) ,1:2:(2*N_z-1))=1;
  %
  V3=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V3((2*N1_x+3):2:(2*N_x+1),         3:2:(2*N1_y+1),1:2:(2*N_z-1))=1;
  V3(         3:2:(2*N_x+1),(2*N1_y+3):2:(2*N_y+1) ,1:2:(2*N_z-1))=1;
  %
  V4=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V4((2*N1_x+1):2:(2*N_x-1),         3:2:(2*N1_y+1),1:2:(2*N_z-1))=1;
  V4(         1:2:(2*N_x-1),(2*N1_y+3):2:(2*N_y+1) ,1:2:(2*N_z-1))=1;
  %
  V5=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V5((2*N1_x+1):2:(2*N_x-1),         1:2:(2*N1_y-1),3:2:(2*N_z+1))=1;
  V5(         1:2:(2*N_x-1),(2*N1_y+1):2:(2*N_y-1) ,3:2:(2*N_z+1))=1;
  %
  V6=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V6((2*N1_x+3):2:(2*N_x+1),         1:2:(2*N1_y-1),3:2:(2*N_z+1))=1;
  V6(         3:2:(2*N_x+1),(2*N1_y+1):2:(2*N_y-1) ,3:2:(2*N_z+1))=1;
  %
  V7=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V7((2*N1_x+3):2:(2*N_x+1),         3:2:(2*N1_y+1),3:2:(2*N_z+1))=1;
  V7(         3:2:(2*N_x+1),(2*N1_y+3):2:(2*N_y+1) ,3:2:(2*N_z+1))=1;
  %
  V8=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V8((2*N1_x+1):2:(2*N_x-1),         3:2:(2*N1_y+1),3:2:(2*N_z+1))=1;
  V8(         1:2:(2*N_x-1),(2*N1_y+3):2:(2*N_y+1) ,3:2:(2*N_z+1))=1;
  
  % logical arrays for midpoints, e.g. V12 represents the midpoints between
  % V1 and V2
  V12=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V12((2*N1_x+2):2:(2*N_x),         1:2:(2*N1_y-1),1:2:(2*N_z-1))=1;
  V12(         2:2:(2*N_x),(2*N1_y+1):2:(2*N_y-1) ,1:2:(2*N_z-1))=1;
  %
  V14=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V14((2*N1_x+1):2:(2*N_x-1),         2:2:(2*N1_y),1:2:(2*N_z-1))=1;
  V14(         1:2:(2*N_x-1),(2*N1_y+2):2:(2*N_y) ,1:2:(2*N_z-1))=1;
  %
  V15=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V15((2*N1_x+1):2:(2*N_x-1),         1:2:(2*N1_y-1),2:2:(2*N_z))=1;
  V15(         1:2:(2*N_x-1),(2*N1_y+1):2:(2*N_y-1) ,2:2:(2*N_z))=1;
  %
  V16=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V16((2*N1_x+2):2:(2*N_x),         1:2:(2*N1_y-1),2:2:(2*N_z))=1;
  V16(         2:2:(2*N_x),(2*N1_y+1):2:(2*N_y-1) ,2:2:(2*N_z))=1;
  %
  V23=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V23((2*N1_x+3):2:(2*N_x+1),         2:2:(2*N1_y),1:2:(2*N_z-1))=1;
  V23(         3:2:(2*N_x+1),(2*N1_y+2):2:(2*N_y) ,1:2:(2*N_z-1))=1;
  %
  V24=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V24((2*N1_x+2):2:(2*N_x),         2:2:(2*N1_y),1:2:(2*N_z-1))=1;
  V24(         2:2:(2*N_x),(2*N1_y+2):2:(2*N_y) ,1:2:(2*N_z-1))=1;
  %
  V26=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V26((2*N1_x+3):2:(2*N_x+1),         1:2:(2*N1_y-1),2:2:(2*N_z))=1;
  V26(         3:2:(2*N_x+1),(2*N1_y+1):2:(2*N_y-1) ,2:2:(2*N_z))=1;
  %
  V34=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V34((2*N1_x+2):2:(2*N_x),         3:2:(2*N1_y+1),1:2:(2*N_z-1))=1;
  V34(         2:2:(2*N_x),(2*N1_y+3):2:(2*N_y+1) ,1:2:(2*N_z-1))=1;
  %
  V36=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V36((2*N1_x+3):2:(2*N_x+1),         2:2:(2*N1_y),2:2:(2*N_z))=1;
  V36(         3:2:(2*N_x+1),(2*N1_y+2):2:(2*N_y) ,2:2:(2*N_z))=1;
  %
  V37=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V37((2*N1_x+3):2:(2*N_x+1),         3:2:(2*N1_y+1),2:2:(2*N_z))=1;
  V37(         3:2:(2*N_x+1),(2*N1_y+3):2:(2*N_y+1) ,2:2:(2*N_z))=1;
  %
  V45=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V45((2*N1_x+1):2:(2*N_x-1),         2:2:(2*N1_y),2:2:(2*N_z))=1;
  V45(         1:2:(2*N_x-1),(2*N1_y+2):2:(2*N_y) ,2:2:(2*N_z))=1;
  %
  V46=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V46((2*N1_x+2):2:(2*N_x),         2:2:(2*N1_y),2:2:(2*N_z))=1;
  V46(         2:2:(2*N_x),(2*N1_y+2):2:(2*N_y) ,2:2:(2*N_z))=1;
  %
  V47=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V47((2*N1_x+2):2:(2*N_x),         3:2:(2*N1_y+1),2:2:(2*N_z))=1;
  V47(         2:2:(2*N_x),(2*N1_y+3):2:(2*N_y+1) ,2:2:(2*N_z))=1;
  %
  V48=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V48((2*N1_x+1):2:(2*N_x-1),         3:2:(2*N1_y+1),2:2:(2*N_z))=1;
  V48(         1:2:(2*N_x-1),(2*N1_y+3):2:(2*N_y+1) ,2:2:(2*N_z))=1;
  %
  V56=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V56((2*N1_x+2):2:(2*N_x),         1:2:(2*N1_y-1),3:2:(2*N_z+1))=1;
  V56(         2:2:(2*N_x),(2*N1_y+1):2:(2*N_y-1) ,3:2:(2*N_z+1))=1;
  %
  V58=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V58((2*N1_x+1):2:(2*N_x-1),        2:2:(2*N1_y),3:2:(2*N_z+1))=1;
  V58(        1:2:(2*N_x-1),(2*N1_y+2):2:(2*N_y) ,3:2:(2*N_z+1))=1;
  %
  V67=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V67((2*N1_x+3):2:(2*N_x+1),         2:2:(2*N1_y),3:2:(2*N_z+1))=1;
  V67(         3:2:(2*N_x+1),(2*N1_y+2):2:(2*N_y) ,3:2:(2*N_z+1))=1;
  %
  V68=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V68((2*N1_x+2):2:(2*N_x),         2:2:(2*N1_y),3:2:(2*N_z+1))=1;
  V68(         2:2:(2*N_x),(2*N1_y+2):2:(2*N_y) ,3:2:(2*N_z+1))=1;
  %
  V78=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V78((2*N1_x+2):2:(2*N_x),         3:2:(2*N1_y+1),3:2:(2*N_z+1))=1;
  V78(         2:2:(2*N_x),(2*N1_y+3):2:(2*N_y+1) ,3:2:(2*N_z+1))=1;

  % used division of the unit cube into 6 tetrahedrons:   
  %   V1 V2 V4 V6 V12 V24 V14 V26 V46 V16
  %   V1 V4 V5 V6 V14 V45 V15 V46 V56 V16
  %   V4 V5 V6 V8 V45 V56 V46 V58 V68 V48
  %   V2 V3 V4 V6 V23 V34 V24 V36 V46 V26
  %   V3 V6 V7 V4 V36 V67 V37 V46 V47 V34
  %   V4 V6 V7 V8 V46 V67 V47 V68 V78 V48
  % size(aux_elem)=(6*10,n_e/6)
  aux_elem=[C(V1)'; C(V2)'; C(V4)'; C(V6)'; C(V12)'; C(V24)'; C(V14)'; C(V26)'; C(V46)'; C(V16)';
            C(V1)'; C(V4)'; C(V5)'; C(V6)'; C(V14)'; C(V45)'; C(V15)'; C(V46)'; C(V56)'; C(V16)';
            C(V4)'; C(V5)'; C(V6)'; C(V8)'; C(V45)'; C(V56)'; C(V46)'; C(V58)'; C(V68)'; C(V48)';
            C(V2)'; C(V3)'; C(V4)'; C(V6)'; C(V23)'; C(V34)'; C(V24)'; C(V36)'; C(V46)'; C(V26)';
            C(V3)'; C(V6)'; C(V7)'; C(V4)'; C(V36)'; C(V67)'; C(V37)'; C(V46)'; C(V47)'; C(V34)';
            C(V4)'; C(V6)'; C(V7)'; C(V8)'; C(V46)'; C(V67)'; C(V47)'; C(V68)'; C(V78)'; C(V48)' ];
        
  % the 10 x n_e array elem      
  elem=reshape(aux_elem,10,n_e);     
  
%
% Surface of the body - the array "surf"
%
  
  % For each face of the body, we define the restriction C_s of the array C 
  % and logical 2D arrays V1_s,...,V24_s which enable to select appropriate
  % nodes from the array C_s. We consider the following ordering of the 
  % nodes within the unit square:
  %      V1_s -> [0 0], V2_s -> [1 0], V3_s -> [1 1], V4_s -> [0 1]
  % V12_s -> [1/2 0], V23_s -> [1 1/2], V34_s -> [1/2 1], V14_s -> [0 1/2], V13_s=V24_s -> [1/2 1/2]
  % Finally, we use the division of a rectangle into 2 triangles which is
  % in accordance to the division of a prism into 6 tetrahedrons, see above.

  % Face 1: z=0 (the bottom of the body)
  C_s=zeros(2*N_x+1,2*N_y+1);
  C_s(:,:)=C(:,:,1);  
   V1_s=false(2*N_x+1,2*N_y+1);  V1_s((2*N1_x+1):2:(2*N_x-1),1:2:(2*N1_y-1))=1;  V1_s(1:2:(2*N_x-1),(2*N1_y+1):2:(2*N_y-1))=1;
   V2_s=false(2*N_x+1,2*N_y+1);  V2_s((2*N1_x+3):2:(2*N_x+1),1:2:(2*N1_y-1))=1;  V2_s(3:2:(2*N_x+1),(2*N1_y+1):2:(2*N_y-1))=1;
   V3_s=false(2*N_x+1,2*N_y+1);  V3_s((2*N1_x+3):2:(2*N_x+1),3:2:(2*N1_y+1))=1;  V3_s(3:2:(2*N_x+1),(2*N1_y+3):2:(2*N_y+1))=1;
   V4_s=false(2*N_x+1,2*N_y+1);  V4_s((2*N1_x+1):2:(2*N_x-1),3:2:(2*N1_y+1))=1;  V4_s(1:2:(2*N_x-1),(2*N1_y+3):2:(2*N_y+1))=1;
  V12_s=false(2*N_x+1,2*N_y+1); V12_s((2*N1_x+2):2:(2*N_x  ),1:2:(2*N1_y-1))=1; V12_s(2:2:(2*N_x  ),(2*N1_y+1):2:(2*N_y-1))=1;
  V23_s=false(2*N_x+1,2*N_y+1); V23_s((2*N1_x+3):2:(2*N_x+1),2:2:(2*N1_y  ))=1; V23_s(3:2:(2*N_x+1),(2*N1_y+2):2:(2*N_y  ))=1;
  V34_s=false(2*N_x+1,2*N_y+1); V34_s((2*N1_x+2):2:(2*N_x  ),3:2:(2*N1_y+1))=1; V34_s(2:2:(2*N_x  ),(2*N1_y+3):2:(2*N_y+1))=1;
  V14_s=false(2*N_x+1,2*N_y+1); V14_s((2*N1_x+1):2:(2*N_x-1),2:2:(2*N1_y  ))=1; V14_s(1:2:(2*N_x-1),(2*N1_y+2):2:(2*N_y  ))=1;
  V24_s=false(2*N_x+1,2*N_y+1); V24_s((2*N1_x+2):2:(2*N_x  ),2:2:(2*N1_y  ))=1; V24_s(2:2:(2*N_x  ),(2*N1_y+2):2:(2*N_y  ))=1;
  aux_surf=[C_s(V1_s)'; C_s(V2_s)'; C_s(V4_s)'; C_s(V24_s)'; C_s(V14_s)'; C_s(V12_s)';
            C_s(V3_s)'; C_s(V4_s)'; C_s(V2_s)'; C_s(V24_s)'; C_s(V23_s)'; C_s(V34_s)'; ];
  surf1=reshape(aux_surf,6,2*(N_x*N_y-N1_x*N1_y));         
  
  % Face 2: z=size_z (the top of the body)
  C_s=zeros(2*N_x+1,2*N_y+1);
  C_s(:,:)=C(:,:,end);  
  aux_surf=[C_s(V1_s)'; C_s(V2_s)'; C_s(V4_s)'; C_s(V24_s)'; C_s(V14_s)'; C_s(V12_s)';
            C_s(V3_s)'; C_s(V4_s)'; C_s(V2_s)'; C_s(V24_s)'; C_s(V23_s)'; C_s(V34_s)'; ];
  surf2=reshape(aux_surf,6,2*(N_x*N_y-N1_x*N1_y));     
  
  % Face 3: y=0 (the front of the body)
  C_s=zeros(2*N_x+1,2*N_z+1);
  C_s(:,:)=C(:,1,:);
   V1_s=false(2*N_x+1,2*N_z+1);   V1_s((2*N1_x+1):2:(2*N_x-1),1:2:(2*N_z-1))=1;
   V2_s=false(2*N_x+1,2*N_z+1);   V2_s((2*N1_x+3):2:(2*N_x+1),1:2:(2*N_z-1))=1;
   V3_s=false(2*N_x+1,2*N_z+1);   V3_s((2*N1_x+3):2:(2*N_x+1),3:2:(2*N_z+1))=1;
   V4_s=false(2*N_x+1,2*N_z+1);   V4_s((2*N1_x+1):2:(2*N_x-1),3:2:(2*N_z+1))=1;
  V12_s=false(2*N_x+1,2*N_z+1);  V12_s((2*N1_x+2):2:(2*N_x  ),1:2:(2*N_z-1))=1;
  V14_s=false(2*N_x+1,2*N_z+1);  V14_s((2*N1_x+1):2:(2*N_x-1),2:2:(2*N_z  ))=1;
  V23_s=false(2*N_x+1,2*N_z+1);  V23_s((2*N1_x+3):2:(2*N_x+1),2:2:(2*N_z  ))=1;
  V13_s=false(2*N_x+1,2*N_z+1);  V13_s((2*N1_x+2):2:(2*N_x  ),2:2:(2*N_z  ))=1;
  V34_s=false(2*N_x+1,2*N_z+1);  V34_s((2*N1_x+2):2:(2*N_x  ),3:2:(2*N_z+1))=1;
  aux_surf=[C_s(V4_s)'; C_s(V1_s)'; C_s(V3_s)'; C_s(V13_s)'; C_s(V34_s)'; C_s(V14_s)';
            C_s(V2_s)'; C_s(V3_s)'; C_s(V1_s)'; C_s(V13_s)'; C_s(V12_s)'; C_s(V23_s)'; ];
  surf3=reshape(aux_surf,6,2*N2_x*N_z);     
  
  % Face 4: x=size_xy (the right hand side of the body)
  C_s=zeros(2*N_y+1,2*N_z+1);
  C_s(:,:)=C(end,:,:);  
   V1_s=false(2*N_y+1,2*N_z+1);   V1_s(1:2:(2*N_y-1),1:2:(2*N_z-1))=1;
   V2_s=false(2*N_y+1,2*N_z+1);   V2_s(3:2:(2*N_y+1),1:2:(2*N_z-1))=1;
   V3_s=false(2*N_y+1,2*N_z+1);   V3_s(3:2:(2*N_y+1),3:2:(2*N_z+1))=1;
   V4_s=false(2*N_y+1,2*N_z+1);   V4_s(1:2:(2*N_y-1),3:2:(2*N_z+1))=1;
  V12_s=false(2*N_y+1,2*N_z+1);  V12_s(2:2:(2*N_y  ),1:2:(2*N_z-1))=1;
  V14_s=false(2*N_y+1,2*N_z+1);  V14_s(1:2:(2*N_y-1),2:2:(2*N_z  ))=1;
  V23_s=false(2*N_y+1,2*N_z+1);  V23_s(3:2:(2*N_y+1),2:2:(2*N_z  ))=1;
  V24_s=false(2*N_y+1,2*N_z+1);  V24_s(2:2:(2*N_y  ),2:2:(2*N_z  ))=1;
  V34_s=false(2*N_y+1,2*N_z+1);  V34_s(2:2:(2*N_y  ),3:2:(2*N_z+1))=1;
  aux_surf=[C_s(V1_s)'; C_s(V2_s)'; C_s(V4_s)'; C_s(V24_s)'; C_s(V14_s)'; C_s(V12_s)';
            C_s(V3_s)'; C_s(V4_s)'; C_s(V2_s)'; C_s(V24_s)'; C_s(V23_s)'; C_s(V34_s)'; ];
  surf4=reshape(aux_surf,6,2*N_y*N_z);       
  
  % Face 5: y=size_xy (the back of the body)
  C_s=zeros(2*N_x+1,2*N_z+1);
  C_s(:,:)=C(:,end,:);
   V1_s=false(2*N_x+1,2*N_z+1);   V1_s(1:2:(2*N_x-1),1:2:(2*N_z-1))=1;
   V2_s=false(2*N_x+1,2*N_z+1);   V2_s(3:2:(2*N_x+1),1:2:(2*N_z-1))=1;
   V3_s=false(2*N_x+1,2*N_z+1);   V3_s(3:2:(2*N_x+1),3:2:(2*N_z+1))=1;
   V4_s=false(2*N_x+1,2*N_z+1);   V4_s(1:2:(2*N_x-1),3:2:(2*N_z+1))=1;
  V12_s=false(2*N_x+1,2*N_z+1);  V12_s(2:2:(2*N_x  ),1:2:(2*N_z-1))=1;
  V14_s=false(2*N_x+1,2*N_z+1);  V14_s(1:2:(2*N_x-1),2:2:(2*N_z  ))=1;
  V23_s=false(2*N_x+1,2*N_z+1);  V23_s(3:2:(2*N_x+1),2:2:(2*N_z  ))=1;
  V13_s=false(2*N_x+1,2*N_z+1);  V13_s(2:2:(2*N_x  ),2:2:(2*N_z  ))=1;
  V34_s=false(2*N_x+1,2*N_z+1);  V34_s(2:2:(2*N_x  ),3:2:(2*N_z+1))=1;
  aux_surf=[C_s(V4_s)'; C_s(V1_s)'; C_s(V3_s)'; C_s(V13_s)'; C_s(V34_s)'; C_s(V14_s)';
            C_s(V2_s)'; C_s(V3_s)'; C_s(V1_s)'; C_s(V13_s)'; C_s(V12_s)'; C_s(V23_s)'; ];
  surf5=reshape(aux_surf,6,2*N_x*N_z);     
  
  % Face 6: x=0 (the left hand side of the body)
  C_s=zeros(2*N_y+1,2*N_z+1);
  C_s(:,:)=C(1,:,:);  
   V1_s=false(2*N_y+1,2*N_z+1);   V1_s((2*N1_y+1):2:(2*N_y-1),1:2:(2*N_z-1))=1;
   V2_s=false(2*N_y+1,2*N_z+1);   V2_s((2*N1_y+3):2:(2*N_y+1),1:2:(2*N_z-1))=1;
   V3_s=false(2*N_y+1,2*N_z+1);   V3_s((2*N1_y+3):2:(2*N_y+1),3:2:(2*N_z+1))=1;
   V4_s=false(2*N_y+1,2*N_z+1);   V4_s((2*N1_y+1):2:(2*N_y-1),3:2:(2*N_z+1))=1;
  V12_s=false(2*N_y+1,2*N_z+1);  V12_s((2*N1_y+2):2:(2*N_y  ),1:2:(2*N_z-1))=1;
  V14_s=false(2*N_y+1,2*N_z+1);  V14_s((2*N1_y+1):2:(2*N_y-1),2:2:(2*N_z  ))=1;
  V23_s=false(2*N_y+1,2*N_z+1);  V23_s((2*N1_y+3):2:(2*N_y+1),2:2:(2*N_z  ))=1;
  V24_s=false(2*N_y+1,2*N_z+1);  V24_s((2*N1_y+2):2:(2*N_y  ),2:2:(2*N_z  ))=1;
  V34_s=false(2*N_y+1,2*N_z+1);  V34_s((2*N1_y+2):2:(2*N_y  ),3:2:(2*N_z+1))=1;
  aux_surf=[C_s(V1_s)'; C_s(V2_s)'; C_s(V4_s)'; C_s(V24_s)'; C_s(V14_s)'; C_s(V12_s)';
            C_s(V3_s)'; C_s(V4_s)'; C_s(V2_s)'; C_s(V24_s)'; C_s(V23_s)'; C_s(V34_s)'; ];
  surf6=reshape(aux_surf,6,2*N2_y*N_z);       
  
  % Face 7: y=size_hole (the hole face perpendicular to y-axis)
  C_s=zeros(2*N_x+1,2*N_z+1);
  C_s(:,:)=C(:,2*N1_y+1,:);
   V1_s=false(2*N_x+1,2*N_z+1);   V1_s(1:2:(2*N1_x-1),1:2:(2*N_z-1))=1;
   V2_s=false(2*N_x+1,2*N_z+1);   V2_s(3:2:(2*N1_x+1),1:2:(2*N_z-1))=1;
   V3_s=false(2*N_x+1,2*N_z+1);   V3_s(3:2:(2*N1_x+1),3:2:(2*N_z+1))=1;
   V4_s=false(2*N_x+1,2*N_z+1);   V4_s(1:2:(2*N1_x-1),3:2:(2*N_z+1))=1;
  V12_s=false(2*N_x+1,2*N_z+1);  V12_s(2:2:(2*N1_x  ),1:2:(2*N_z-1))=1;
  V14_s=false(2*N_x+1,2*N_z+1);  V14_s(1:2:(2*N1_x-1),2:2:(2*N_z  ))=1;
  V23_s=false(2*N_x+1,2*N_z+1);  V23_s(3:2:(2*N1_x+1),2:2:(2*N_z  ))=1;
  V13_s=false(2*N_x+1,2*N_z+1);  V13_s(2:2:(2*N1_x  ),2:2:(2*N_z  ))=1;
  V34_s=false(2*N_x+1,2*N_z+1);  V34_s(2:2:(2*N1_x  ),3:2:(2*N_z+1))=1;
  aux_surf=[C_s(V4_s)'; C_s(V1_s)'; C_s(V3_s)'; C_s(V13_s)'; C_s(V34_s)'; C_s(V14_s)';
            C_s(V2_s)'; C_s(V3_s)'; C_s(V1_s)'; C_s(V13_s)'; C_s(V12_s)'; C_s(V23_s)'; ];
  surf7=reshape(aux_surf,6,2*N1_x*N_z);            
  
  % Face 8: x=size_hole (the hole face perpendicular to x-axis)
  C_s=zeros(2*N_y+1,2*N_z+1);
  C_s(:,:)=C(2*N1_x+1,:,:);  
   V1_s=false(2*N_y+1,2*N_z+1);   V1_s(1:2:(2*N1_y-1),1:2:(2*N_z-1))=1;
   V2_s=false(2*N_y+1,2*N_z+1);   V2_s(3:2:(2*N1_y+1),1:2:(2*N_z-1))=1;
   V3_s=false(2*N_y+1,2*N_z+1);   V3_s(3:2:(2*N1_y+1),3:2:(2*N_z+1))=1;
   V4_s=false(2*N_y+1,2*N_z+1);   V4_s(1:2:(2*N1_y-1),3:2:(2*N_z+1))=1;
  V12_s=false(2*N_y+1,2*N_z+1);  V12_s(2:2:(2*N1_y  ),1:2:(2*N_z-1))=1;
  V14_s=false(2*N_y+1,2*N_z+1);  V14_s(1:2:(2*N1_y-1),2:2:(2*N_z  ))=1;
  V23_s=false(2*N_y+1,2*N_z+1);  V23_s(3:2:(2*N1_y+1),2:2:(2*N_z  ))=1;
  V24_s=false(2*N_y+1,2*N_z+1);  V24_s(2:2:(2*N1_y  ),2:2:(2*N_z  ))=1;
  V34_s=false(2*N_y+1,2*N_z+1);  V34_s(2:2:(2*N1_y  ),3:2:(2*N_z+1))=1;
  aux_surf=[C_s(V1_s)'; C_s(V2_s)'; C_s(V4_s)'; C_s(V24_s)'; C_s(V14_s)'; C_s(V12_s)';
            C_s(V3_s)'; C_s(V4_s)'; C_s(V2_s)'; C_s(V24_s)'; C_s(V23_s)'; C_s(V34_s)'; ];
  surf8=reshape(aux_surf,6,2*N1_y*N_z);       
  
  % the array "surf"
  surf = [surf1 surf2 surf3 surf4 surf5 surf6 surf7 surf8] ;
  
  
%
% Boundary conditions
%
  
  % Nonhomogeneous Neumann boundary conditions on Face 5
  neumann=surf5;     

  % logical array indicating the nodes with the Dirichlet boundary cond.
  Q = coord>0 ;
  
end
