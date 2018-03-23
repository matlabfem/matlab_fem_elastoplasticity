function [dshape] = shapeder (point,etype)
% SHAPEDER Returns the gradients of the shape functions with
%          respect to the reference coordinates (xi,eta,...).
%
%  point : point(nod,nop), the coordinates of the
%          points on the reference element.
% dshape : dshape(nod,nos,nop), the gradients of the
%          shape functions (second) at all points (third)
%          with respect to the reference cordinates.
%  etype : 'P0','P1','P2', etc., the element type.
%         
%          Note: 
%          nod - dimension of the element.
%          nop - number of points.
%          nos - number of shape functions.
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
 
          dshape = zeros(1,1,nop);

     case {'P1'}
     % Linear shape functions.
  
          dshape = [1 -1];
          dshape = reshape(dshape,2,1);
          dshape = dshape*ones(1,nop);
          dshape = reshape(dshape,1,2,nop);
  
     case {'P2'}
  
          dshape = [  4 .* l1 - 1;    ...
                    - 4 .* l2 + 1;    ...
                      4 .* (l2 - l1)];
          dshape = reshape(dshape,1,3,nop);

     otherwise, error('Only P1, P2 elements implemented.');
  
     end
  
case {2}
% 2-D elements.

     l1 = point(1,:);
     l2 = point(2,:);
     l3 = 1 - l1 - l2;
   
     switch etype
   
     case {'P0'}
     % Constant shape function.

          dshape = zeros(2,1,nop);

     case {'P1'},
     % Linear shape functions.
   
          dshape = [1 0 -1;
                    0 1 -1];
          dshape = reshape(dshape,6,1);
          dshape = dshape*ones(1,nop);
          dshape = reshape(dshape,2,3,nop);
   
     case {'P2'}
         %comment by Jan Valdman
         %l1=x, l2=y, l3=1-x-y are barycentric coordinates
         %6 basis functions
         %on point (0,0)        phi1=(2*l3-1)*l3   node 0
         %on point (1,0)        phi2=(2*l1-1)*l1   node 1
         %on point (0,1)        phi3=(2*l2-1)*l1   node 2
         %on point (0.5,0.5)    phi4=4*l1*l2       edge 12
         %on point (0,0.5)      phi5=4*l2*l3       edge 23
         %on point (0.5,0)      phi6=4*l1*l3       edge 13
         %note that d(l3)/d(l1)=-l1, d(l3)/d(l2)=-l2

          dshape = [  4 .* l1 - 1;    ...       %d(phi2)/dx
                      zeros(1,nop);   ...       %d(phi2)/dy
                      zeros(1,nop);   ...       %d(phi3)/dx
                      4 .* l2 - 1;    ...       %d(phi3)/dy
                    - 4 .* l3 + 1;    ...       %d(phi1)/dx
                    - 4 .* l3 + 1;    ...       %d(phi1)/dy
                      4 .* l2;        ...       %d(phi4)/dx
                      4 .* l1;        ...       %d(phi4)/dy
                    - 4 .* l2;        ...       %d(phi5)/dx 
                      4 .* (l3 - l2); ...       %d(phi5)/dy
                      4 .* (l3 - l1); ...       %d(phi6)/dx
                    - 4 .* l1];                 %d(phi6)/dy
                
           dshape = [ - 4 .* l3 + 1;    ...     %d(phi1)/dx
                    - 4 .* l3 + 1;    ...       %d(phi1)/dy
                      4 .* l1 - 1;    ...       %d(phi2)/dx
                      zeros(1,nop);   ...       %d(phi2)/dy
                      zeros(1,nop);   ...       %d(phi3)/dx
                      4 .* l2 - 1;    ...       %d(phi3)/dy                
                      4 .* l2;        ...       %d(phi4)/dx
                      4 .* l1;        ...       %d(phi4)/dy
                    - 4 .* l2;        ...       %d(phi5)/dx 
                      4 .* (l3 - l2); ...       %d(phi5)/dy
                      4 .* (l3 - l1); ...       %d(phi6)/dx
                    - 4 .* l1];                 %d(phi6)/dy     
                 
          dshape = reshape(dshape,2,6,nop);

     otherwise, error('Only P1 and P2 elements implemented.');
   
     end

case {3}
% 3-D elements.

     l1 = point(1,:);
     l2 = point(2,:);
     l3 = point(3,:);
     l4 = 1 - l1 - l2 -l3;
   
     switch etype
   
     case {'P0'}
     % Constant shape function.
 
          dshape = zeros(1,1,nop);

     case {'P1'}
     % Linear shape functions.
   
          dshape = [1 0 0 -1;
                    0 1 0 -1;
                    0 0 1 -1];
          dshape = reshape(dshape,12,1);
          dshape = dshape*ones(1,nop);
          dshape = reshape(dshape,3,4,nop);
     
     case {'P2'}
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
     
     
     dshape =       [ - 4 .* l4 + 1;    ...     %d(phi1)/dx
                      - 4 .* l4 + 1;    ...     %d(phi1)/dy
                      - 4 .* l4 + 1;    ...     %d(phi1)/dz
                        4 .* l1 - 1;    ...     %d(phi2)/dx
                        zeros(1,nop);   ...     %d(phi2)/dy
                        zeros(1,nop);   ...     %d(phi2)/dz        
                        zeros(1,nop);   ...     %d(phi3)/dx
                        4 .* l2 - 1;    ...     %d(phi3)/dy  
                        zeros(1,nop);   ...     %d(phi3)/dz 
                        zeros(1,nop);   ...     %d(phi4)/dx
                        zeros(1,nop);   ...     %d(phi4)/dy 
                        4 .* l3 - 1;    ...     %d(phi4)/dz  
                        4 .* (l4 - l1); ...     %d(phi5)/dx
                      - 4 .* l1;        ...     %d(phi5)/dy
                      - 4 .* l1;        ...     %d(phi5)/dz
                      - 4 .* l2;        ...     %d(phi6)/dx
                        4 .* (l4 - l2); ...     %d(phi6)/dy
                      - 4 .* l2;        ...     %d(phi6)/dz
                      - 4 .* l3;        ...     %d(phi7)/dx               
                      - 4 .* l3;        ...     %d(phi7)/dy
                        4 .* (l4 - l3); ...     %d(phi7)/dz
                        4 .* l2;        ...     %d(phi8)/dx
                        4 .* l1;        ...     %d(phi8)/dy
                        zeros(1,nop);   ...     %d(phi8)/dz 
                        4 .* l3;        ...     %d(phi9)/dx
                        zeros(1,nop);   ...     %d(phi9)/dy 
                        4 .* l1;        ...     %d(phi9)/dz
                        zeros(1,nop);   ...     %d(phi10)/dx 
                        4 .* l3;        ...     %d(phi10)/dx
                        4 .* l2];               %d(phi10)/dz
                                   
      dshape = reshape(dshape,3,10,nop);

     case {'Q1'}
     % Tri-linear (8 node) hexahedron
     % (x,y,z) is a point on local coordinate
     % 8 basis nodal functions
     x = point(1,:);
     y = point(2,:);
     z = point(3,:);

     % Look at page 3, paper "Simulation of Hyperelastic Materials Using Energy Constraints"
     
     s = [0 0 0;                % ATTENTION ON ORDER OF POINTS IN S
          1 0 0;
          1 1 0;
          0 1 0;
          0 0 1; 
          1 0 1;
          1 1 1;
          0 1 1] * 2 - 1;

     dshape = zeros(3,8,nop);   % 3: 3D x,y,z; 8: basis of 8 functions; nop: number of points
     
     for i = 1:8
         dshape(1,i,:) = 1/8 * s(i,1) .* (1+s(i,2)*y) .* (1+s(i,3)*z);  % d(phi_i)/dx
         dshape(2,i,:) = 1/8 * (1+s(i,1)*x) .* s(i,2) .* (1+s(i,3)*z);  % d(phi_i)/dy
         dshape(3,i,:) = 1/8 * (1+s(i,1)*x) .* (1+s(i,2)*y) .* s(i,3);  % d(phi_i)/dz
     end
     
     otherwise, error('Only P1, P2, and Q1 elements implemented.');
     
     end

end

return
