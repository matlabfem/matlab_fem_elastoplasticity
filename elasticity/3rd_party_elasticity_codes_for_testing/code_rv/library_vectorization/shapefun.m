function [shape] = shapefun (point,etype)
% SHAPEFUN Returns the values of the shape functions.
%
% point : point(nod,nop), the coordinates of the
%         points on the reference element.
% shape : shape(nos,nop), the values of the shape 
%         functions (first) at all points (second).
% etype : 'P0','P1','P2', etc., the element type.
%         
%         Note: 
%         nod - dimension of the elements.
%         nop - number of points.
%         nos - number of shape functions.
%
% Quadratic Q1 added

nod = size(point,1);
nop = size(point,2);

switch nod

case {1},
% 1-D elements.

     l1 = point(1,:);
     l2 = 1-l1;
  
     switch etype
  
     case {'P0'}
     % Constant shape function.

          shape = ones(1,nop);
          
     case {'P1'},
     % Linear shape functions.
  
          shape = [l1; ...
                   l2];

     case {'P2'},
     % Quadratic shape functions.

       shape = [(2 .* l1 - 1) .* l1; ...
                (2 .* l2 - 1) .* l2; ...
                4 .* l1 .* l2];

     otherwise, error('Only P1, P2 elements implemented.');

     end

case {2},
% 2-D elements.

     l1 = point(1,:);
     l2 = point(2,:);
     l3 = 1 - l1 - l2;

     switch etype

     case {'P0'}
     % Constant shape function.

          shape = ones(1,nop);

     case {'P1'},
     % Linear shape functions.

          shape = [l1; ...
                   l2; ...
                   l3];
    
     case {'P2'},
     % Quadratic shape functions.

          shape = [(2 .* l1 - 1) .* l1; ...
                   (2 .* l2 - 1) .* l2; ...
                   (2 .* l3 - 1) .* l3; ... 
                   4 .* l1 .* l2;       ...
                   4 .* l2 .* l3;       ...
                   4 .* l3 .* l1];
               
          %modified by Jan Valdman     
          shape = [(2 .* l3 - 1) .* l3; ... 
                   (2 .* l1 - 1) .* l1; ...
                   (2 .* l2 - 1) .* l2; ...                
                   4 .* l1 .* l2;       ...
                   4 .* l2 .* l3;       ...
                   4 .* l3 .* l1];     
               

     otherwise, error('Only P1, P2 elements implemented.');

     end

case {3}, 
% 3-D elements.

     l1 = point(1,:);
     l2 = point(2,:);
     l3 = point(3,:);
     l4 = 1 - l1 - l2 -l3;

     switch etype
 
     case {'P0'}
     % Constant shape function.

          shape = ones(1,nop);

     case {'P1'},
     % Linear shape functions.

          shape = [l1; ...
                   l2; ...
                   l3; ...
                   l4];
               
      case {'P2'},

     %Quadratic shape functions
     %l1=x, l2=y, l3=z, l4=1-x-y-z are barycentric coordinates
     %10 basis nodal functions
     %on point (0,0,0)        phi1=(2*l4-1)*l4         node 0
     %on point (1,0,0)        phi2=(2*l1-1)*l1         node 1
     %on point (0,1,0)        phi3=(2*l2-1)*l2         node 2      
     %on point (0,0,1)        phi4=(2*l3-1)*l3         node 3
     %on point (0.5,0,0)      phi5=4*l1*l4             edge 01
     %on point (0,0.5,0)      phi6=4*l2*l4             edge 02
     %on point (0,0,0.5)      phi7=4*l3*l4             edge 03
     %on point (0.5,0.5,0)    phi8=4*l1*l2             edge 12
     %on point (0.5,0,0.5)    phi9=4*l1*l3             edge 13
     %on point (0,0.5,0.5)    phi10=4*l2*l3            edge 23
     
     
     shape =       [  (2 .* l4 - 1) .* l4; ...
                      (2 .* l1 - 1) .* l1; ...
                      (2 .* l2 - 1) .* l2; ...
                      (2 .* l3 - 1) .* l3; ...
                       4 .* l1 .* l4;       ...
                       4 .* l2 .* l4;       ...
                       4 .* l3 .* l4;       ...
                       4 .* l1 .* l2;       ...
                       4 .* l1 .* l3;       ...
                       4 .* l2 .* l3; ] ;        
          
     
     case {'Q1'},
     % Tri-linear (8 node) hexahedron
     % (x,y,z) is a point on local coordinate
     % 8 basis nodal functions
     x = point(1,:);
     y = point(2,:);
     z = point(3,:);

     % Look at page 3, paper "Simulation of Hyperelastic Materials Using Energy Constraints"
     shape = zeros(8,nop);

     s = [0 0 0;                % ATTENTION ON ORDER OF POINTS IN S
          1 0 0;
          1 1 0;
          0 1 0;
          0 0 1; 
          1 0 1;
          1 1 1;
          0 1 1] * 2 - 1;
     
     for i = 1:8
         shape(i,:) = 1/8 * (1+s(i,1)*x) .* (1+s(i,2)*y) .* (1+s(i,3)*z);
     end
 
     otherwise, error('Only P1, P2, and Q1 elements implemented.');

     end

end

return
