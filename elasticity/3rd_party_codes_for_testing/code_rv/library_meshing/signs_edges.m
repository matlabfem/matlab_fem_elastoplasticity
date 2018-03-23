function signs = signs_edges(elems2nodes)

% SIGNS_EDGES
%   Calculate signs for the edges of a 2D/3D mesh.
%
%   This data is needed in order to use linear Nedelec's elements in 2D/3D
%   and the linear Raviart-Thomas element in 2D. For the linear Raviart-
%   Thomas element in 3D we need signs related to faces, and this data is
%   provided by the function 'signs_faces(...)'.
%
%   The edge signs can be easily deduced from the mesh data itself by
%   directly using the data structure which represents the elements
%   (triangles or tetrahedrons) by their nodes. The signs are obtained
%   with minimal matrix operations in a vectorized manner.
%
% SYNTAX:  signs = signs_edges(elements)
%
% IN:   elems2nodes  Elements by their nodes.
%                    In 2D, elements(i,1:3) are the three nodes of the i'th
%                    triangle, and in 3D elements(i,1:6) are the six nodes of
%                    the i'th tetrahedron.
%
% OUT:  signs        Signs for element edges, corresponding to the data
%                    structure 'elements': signs(i,j) is the sign related to
%                    the j'th edge of the i'th triangle (or tetrahedron).
%

dim = size(elems2nodes,2)-1;

if ( dim == 2 )
    tmp = elems2nodes(:,[2 3 1]) - elems2nodes(:,[3 1 2]);
    signs = tmp ./ abs(tmp);
elseif (dim == 3)
    tmp = elems2nodes(:,[1 1 1 2 3 4]) - elems2nodes(:,[2 3 4 3 4 2]);
    signs = tmp ./ abs(tmp);
else
    error('The data is not understood.')
end

end