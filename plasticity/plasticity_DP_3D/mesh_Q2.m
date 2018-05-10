
function [coord,elem,surf,dirichlet,Q]=mesh_Q2(level,size_xy,size_z)

% =========================================================================
%
%  This function creates tetrahedral mesh for Q2 elements consisting of 
%  8 vertices and 12 midpoints.
%
%  input data:
%    level   - an integer defining a density of a uniform mesh
%    size_xy - size of the body in directions x and y (integer)
%    size_z  - size of the body in z-direction (integer) 
%    body=(0,size_xy)x(0,size_xy)x(0,size_z)
%
%  output data:
%    coord     - coordinates of the nodes, size(coord)=(3,n_n) where n_n is a
%                number of nodes including midpoints
%    elem      - array containing numbers of nodes defining each element, 
%                size(elem)=(20,n_e), n_e = number of elements
%    surf      - array containing numbers of nodes defining each surface element, 
%                size(surf)=(8,n_s), n_s = number of surface elements
%    dirichlet - array indicating the nodes where the nonhomogenous
%                Dirichlet boundary condition is considered,
%                size(dirichlet)=(3,n_n)
%    Q         - logical array indicating the nodes where the Dirichlet
%                boundary condition is considered, size(Q)=(3,n_n)
%
% ======================================================================
%

%
% numbers of segments, nodes and elements
%

  N_x = size_xy*2^level;      % number of segments in x direction
  N_y = N_x;                  % number of segments in y direction
  N_z = size_z*2^level;       % number of segments in z-direction
   
%
% C - 3D auxilliary array that contains node numbers and that is important 
% for the mesh construction. To specify selected midpoints, we define 
% a 3D logical array Q_C.
%
  C=zeros(2*N_x+1,2*N_y+1,2*N_z+1);
  %
  Q_C=true(2*N_x+1,2*N_y+1,2*N_z+1);
  Q_C(1:2:(2*N_x-1),2:2:(2*N_y  ),2:2:(2*N_z  ))=0;
  Q_C(2:2:(2*N_x  ),3:2:(2*N_y+1),2:2:(2*N_z  ))=0;
  Q_C(3:2:(2*N_x+1),2:2:(2*N_y  ),2:2:(2*N_z  ))=0;
  Q_C(2:2:(2*N_x  ),1:2:(2*N_y-1),2:2:(2*N_z  ))=0;
  Q_C(2:2:(2*N_x  ),2:2:(2*N_y  ),1:2:(2*N_z-1))=0;
  Q_C(2:2:(2*N_x  ),2:2:(2*N_y  ),3:2:(2*N_z+1))=0;
  Q_C(2:2:(2*N_x  ),2:2:(2*N_y  ),2:2:(2*N_z  ))=0;
  %
  C(Q_C)=1:length(C(Q_C));
  
%
% coordinates of nodes
%
  % coordinates in directions x, y and z
  coord_x=linspace(0,size_xy,2*N_x+1);
  coord_y=linspace(0,size_xy,2*N_y+1);
  coord_z=linspace(0,size_z,2*N_z+1);
  
  % 3D arrays containing coordinates in x,y and z directions
  c_x=repmat(coord_x,[1,(2*N_y+1),(2*N_z+1)]);
  %
  c_y=repmat(coord_y,[2*N_x+1,1,(2*N_z+1)]);
  %
  c_z=reshape(kron(coord_z,ones(1,(2*N_x+1)*(2*N_y+1))),[(2*N_x+1),(2*N_y+1),(2*N_z+1)]);  
  
  % the required array of coordinates, size(coord)=(3,n_n)
  coord=[c_x(Q_C)'; c_y(Q_C)'; c_z(Q_C)'] ;
  
% 
% construction of the array elem, size(elem)=(20,n_e)
%
  % ordering of the nodes creating the unit cube:
  %  V1 -> [0 0 0], V2 -> [1 0 0], V3 -> [1 1 0], V4 -> [0 1 0]
  %  V5 -> [0 0 1], V6 -> [1 0 1], V7 -> [1 1 1], V8 -> [0 1 1]
  %  V1,...,V8 are logical 3D arrays which enable to select appropriate
  %  nodes from the array C.

  V1=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V1(1:2:(2*N_x-1),1:2:(2*N_y-1),1:2:(2*N_z-1))=1;
  %
  V2=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V2(3:2:(2*N_x+1),1:2:(2*N_y-1),1:2:(2*N_z-1))=1;
  %
  V3=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V3(3:2:(2*N_x+1),3:2:(2*N_y+1),1:2:(2*N_z-1))=1;
  %
  V4=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V4(1:2:(2*N_x-1),3:2:(2*N_y+1),1:2:(2*N_z-1))=1;
  %
  V5=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V5(1:2:(2*N_x-1),1:2:(2*N_y-1),3:2:(2*N_z+1))=1;
  %
  V6=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V6(3:2:(2*N_x+1),1:2:(2*N_y-1),3:2:(2*N_z+1))=1;
  %
  V7=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V7(3:2:(2*N_x+1),3:2:(2*N_y+1),3:2:(2*N_z+1))=1;
  %
  V8=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V8(1:2:(2*N_x-1),3:2:(2*N_y+1),3:2:(2*N_z+1))=1;
  
  % logical arrays for midpoints, e.g. V12 represents the midpoints between
  % V1 and V2
  V12=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V12(2:2:(2*N_x),1:2:(2*N_y-1),1:2:(2*N_z-1))=1;
  %
  V14=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V14(1:2:(2*N_x-1),2:2:(2*N_y),1:2:(2*N_z-1))=1;
  %
  V15=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V15(1:2:(2*N_x-1),1:2:(2*N_y-1),2:2:(2*N_z))=1;
  %
  V23=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V23(3:2:(2*N_x+1),2:2:(2*N_y),1:2:(2*N_z-1))=1;
  %
  V26=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V26(3:2:(2*N_x+1),1:2:(2*N_y-1),2:2:(2*N_z))=1;
  %
  V34=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V34(2:2:(2*N_x),3:2:(2*N_y+1),1:2:(2*N_z-1))=1;
  %
  V37=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V37(3:2:(2*N_x+1),3:2:(2*N_y+1),2:2:(2*N_z))=1;
  %
  V48=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V48(1:2:(2*N_x-1),3:2:(2*N_y+1),2:2:(2*N_z))=1;
  %
  V56=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V56(2:2:(2*N_x),1:2:(2*N_y-1),3:2:(2*N_z+1))=1;
  %
  V58=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V58(1:2:(2*N_x-1),2:2:(2*N_y),3:2:(2*N_z+1))=1;
  %
  V67=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V67(3:2:(2*N_x+1),2:2:(2*N_y),3:2:(2*N_z+1))=1;
  %
  V78=false(2*N_x+1,2*N_y+1,2*N_z+1);
  V78(2:2:(2*N_x),3:2:(2*N_y+1),3:2:(2*N_z+1))=1;
  
  % used ordering of nodes within a cube:   
  % V1 V2 V3 V4 V5 V6 V7 V8 V12 V23 V34 V14 V56 V67 V78 V58 V15 V26 V37 V48
  
  elem=[ C(V1)';  C(V2)';  C(V3)';  C(V4)'; 
         C(V5)';  C(V6)';  C(V7)';  C(V8)'; 
        C(V12)'; C(V23)'; C(V34)'; C(V14)'; 
        C(V56)'; C(V67)'; C(V78)'; C(V58)'; 
        C(V15)'; C(V26)'; C(V37)'; C(V48)' ];    
    
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
   V1_s=false(2*N_x+1,2*N_y+1);  V1_s(1:2:(2*N_x-1),1:2:(2*N_y-1))=1;  
   V2_s=false(2*N_x+1,2*N_y+1);  V2_s(3:2:(2*N_x+1),1:2:(2*N_y-1))=1;  
   V3_s=false(2*N_x+1,2*N_y+1);  V3_s(3:2:(2*N_x+1),3:2:(2*N_y+1))=1;  
   V4_s=false(2*N_x+1,2*N_y+1);  V4_s(1:2:(2*N_x-1),3:2:(2*N_y+1))=1;  
  V12_s=false(2*N_x+1,2*N_y+1); V12_s(2:2:(2*N_x  ),1:2:(2*N_y-1))=1; 
  V23_s=false(2*N_x+1,2*N_y+1); V23_s(3:2:(2*N_x+1),2:2:(2*N_y  ))=1; 
  V34_s=false(2*N_x+1,2*N_y+1); V34_s(2:2:(2*N_x  ),3:2:(2*N_y+1))=1; 
  V14_s=false(2*N_x+1,2*N_y+1); V14_s(1:2:(2*N_x-1),2:2:(2*N_y  ))=1; 
  surf1=[ C_s(V1_s)';  C_s(V2_s)';  C_s(V3_s)';  C_s(V4_s)'; 
          C_s(V12_s)'; C_s(V23_s)'; C_s(V34_s)'; C_s(V14_s)' ];
  
  % Face 2: z=size_z (the top of the body)
  C_s=zeros(2*N_x+1,2*N_y+1);
  C_s(:,:)=C(:,:,end);  
  surf2=[ C_s(V1_s)';  C_s(V2_s)';  C_s(V3_s)';  C_s(V4_s)'; 
          C_s(V12_s)'; C_s(V23_s)'; C_s(V34_s)'; C_s(V14_s)' ];
  
  % Face 3: y=0 (the front of the body)
  C_s=zeros(2*N_x+1,2*N_z+1);
  C_s(:,:)=C(:,1,:);
   V1_s=false(2*N_x+1,2*N_z+1);   V1_s(1:2:(2*N_x-1),1:2:(2*N_z-1))=1;
   V2_s=false(2*N_x+1,2*N_z+1);   V2_s(3:2:(2*N_x+1),1:2:(2*N_z-1))=1;
   V3_s=false(2*N_x+1,2*N_z+1);   V3_s(3:2:(2*N_x+1),3:2:(2*N_z+1))=1;
   V4_s=false(2*N_x+1,2*N_z+1);   V4_s(1:2:(2*N_x-1),3:2:(2*N_z+1))=1;
  V12_s=false(2*N_x+1,2*N_z+1);  V12_s(2:2:(2*N_x  ),1:2:(2*N_z-1))=1;
  V14_s=false(2*N_x+1,2*N_z+1);  V14_s(1:2:(2*N_x-1),2:2:(2*N_z  ))=1;
  V23_s=false(2*N_x+1,2*N_z+1);  V23_s(3:2:(2*N_x+1),2:2:(2*N_z  ))=1;
  V34_s=false(2*N_x+1,2*N_z+1);  V34_s(2:2:(2*N_x  ),3:2:(2*N_z+1))=1;
  surf3=[ C_s(V1_s)';  C_s(V2_s)';  C_s(V3_s)';  C_s(V4_s)'; 
          C_s(V12_s)'; C_s(V23_s)'; C_s(V34_s)'; C_s(V14_s)' ];
  
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
  V34_s=false(2*N_y+1,2*N_z+1);  V34_s(2:2:(2*N_y  ),3:2:(2*N_z+1))=1;
  surf4=[ C_s(V1_s)';  C_s(V2_s)';  C_s(V3_s)';  C_s(V4_s)'; 
          C_s(V12_s)'; C_s(V23_s)'; C_s(V34_s)'; C_s(V14_s)' ];
  
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
  V34_s=false(2*N_x+1,2*N_z+1);  V34_s(2:2:(2*N_x  ),3:2:(2*N_z+1))=1;
  surf5=[ C_s(V1_s)';  C_s(V2_s)';  C_s(V3_s)';  C_s(V4_s)'; 
          C_s(V12_s)'; C_s(V23_s)'; C_s(V34_s)'; C_s(V14_s)' ];
  
  % Face 6: x=0 (the left hand side of the body)
  C_s=zeros(2*N_y+1,2*N_z+1);
  C_s(:,:)=C(1,:,:);  
   V1_s=false(2*N_y+1,2*N_z+1);   V1_s(1:2:(2*N_y-1),1:2:(2*N_z-1))=1;
   V2_s=false(2*N_y+1,2*N_z+1);   V2_s(3:2:(2*N_y+1),1:2:(2*N_z-1))=1;
   V3_s=false(2*N_y+1,2*N_z+1);   V3_s(3:2:(2*N_y+1),3:2:(2*N_z+1))=1;
   V4_s=false(2*N_y+1,2*N_z+1);   V4_s(1:2:(2*N_y-1),3:2:(2*N_z+1))=1;
  V12_s=false(2*N_y+1,2*N_z+1);  V12_s(2:2:(2*N_y  ),1:2:(2*N_z-1))=1;
  V14_s=false(2*N_y+1,2*N_z+1);  V14_s(1:2:(2*N_y-1),2:2:(2*N_z  ))=1;
  V23_s=false(2*N_y+1,2*N_z+1);  V23_s(3:2:(2*N_y+1),2:2:(2*N_z  ))=1;
  V34_s=false(2*N_y+1,2*N_z+1);  V34_s(2:2:(2*N_y  ),3:2:(2*N_z+1))=1;
  surf6=[ C_s(V1_s)';  C_s(V2_s)';  C_s(V3_s)';  C_s(V4_s)'; 
          C_s(V12_s)'; C_s(V23_s)'; C_s(V34_s)'; C_s(V14_s)' ];  
   
  % the array "surf"
  surf = [surf1 surf2 surf3 surf4 surf5 surf6] ;    
  
%
% Boundary conditions
%
    
  % array indicating the nodes with non-homogen. Dirichlet boundary cond.
  dirichlet = zeros(size(coord));
  dirichlet(2,(coord(2,:)==size_xy)&(coord(1,:)<=1.0001)) = 1;     
  
  % logical array indicating the nodes with the Dirichlet boundary cond.
  Q = coord>0 ;
  Q(2,(coord(2,:)==size_xy)&(coord(1,:)<=1.0001)) = 0;      
  Q(3,(coord(3,:)==size_z)) = 0;     
  Q(1,(coord(1,:)==size_xy)) = 0;   
   
end
