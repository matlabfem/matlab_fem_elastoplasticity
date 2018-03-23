function bfaces2nodes = get_boundary_faces(elems2faces,faces2nodes)

% GETBOUNDARYFACES
%   Calculates the boundary faces by their nodes for 3D mesh. In 2D this
%   data is calculated by 'refinement_uniform()'.
%
% IN:  elems2faces      elements by their faces
%      face2nodes       faces by their nodes
%
% OUT: bfaces2nodes     boundary faces by their nodes
%

is_quad = 0;
n_col_E = 8;
col_E_ind = 5:8;

if size(elems2faces,2) == 6
  is_quad = 1;
  n_col_E = 12;
  col_E_ind = 7:12;
end
  

[A1,I1] = sort(elems2faces(:,1));
[A2,I2] = sort(elems2faces(:,2));
[A3,I3] = sort(elems2faces(:,3));
[A4,I4] = sort(elems2faces(:,4));
if is_quad
  [A5,I5] = sort(elems2faces(:,5));
  [A6,I6] = sort(elems2faces(:,6));
end

nfaces = max(max(elems2faces));
E = zeros(nfaces,n_col_E);

E(A1,1) = I1;
E(A2,2) = I2;
E(A3,3) = I3;
E(A4,4) = I4;
if is_quad
  E(A5,5) = I5;
  E(A6,6) = I6;  
end

% If the same face is listed in the same row of 'elems2faces' more than,
% once it will simply be missed! Because of this we have to insert the
% following dummy variables in order to determine the boundary faces.
ind1 = (diff(A1) == 0);
ind2 = (diff(A2) == 0);
ind3 = (diff(A3) == 0);
ind4 = (diff(A4) == 0);
if is_quad
  ind5 = (diff(A5) == 0);
  ind6 = (diff(A6) == 0);  
end

E(A1(ind1),col_E_ind(1)) = 1;
E(A2(ind2),col_E_ind(2)) = 1;
E(A3(ind3),col_E_ind(3)) = 1;
E(A4(ind4),col_E_ind(4)) = 1;
if is_quad
  E(A5(ind5),col_E_ind(5)) = 1;
  E(A6(ind6),col_E_ind(6)) = 1;  
end

% final sorting
E = sort(E,2,'descend');

% Get boundary nodes by first examining which columns in E
% have only one nonzero element, meaning that this face is
% related to only one single tetra, which means it is on the
% boundary of the domain. Since faces are defined by their nodes,
% we have the boundary nodes too.
ind = (E(:,2) == 0);
bfaces2nodes = faces2nodes(ind,:);

end
