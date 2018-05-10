function signs = signs_faces(nodes2coord,elems2faces,faces2nodes,B_K)

% SIGNS_FACES
%   Calculate signs for the faces of a 3D mesh (tetrahedrons).
%
%   This data is needed for the linear Raviart-Thomas elements in 3D. In
%   2D the signs are related to edges, and the edge signs are calculated
%   by the function 'signs_edges(...)'.
%
%   Unlike sign calculation for edges, the orientation of faces cannot
%   be deduced directly from the mesh data. Instead, we must first choose
%   a positive normal direction for each face. Then we simply calculate
%   the outward normals for each tetrahedron's face, and compare it to the
%   chosen positive normal direction.
%
%   The face signs are obtained with matrix operations in a vectorized
%   manner. However, the number of matrix operations is much more higher
%   than in edge sign calculation, making it considerably slower.
%
% SYNTAX:  signs = signs_faces(coordinates, elems2faces, faces2nodes, B_K)
%
% IN:   nodes2coord    nodes by coordinates
%       elems2faces    elements by their faces
%       faces2nodes    faces by their nodes
%       B_K            the matrix part of the affine transformations
%
% OUT:  signs          signs for element faces, corresponding to the data
%                      structure 'elems2faces': signs(i,j) is the sign
%                      related to the j'th face of the i'th tetrahedron.
%

dim = size(B_K,1);
if ( dim ~= 3 )
    error('Face signs can be calculated only in 3D.')
end

% calculate normal vector for each face, and let us choose this
% normal direction as the positive direction
p1 = nodes2coord(faces2nodes(:,1),1:3);     % first,
p2 = nodes2coord(faces2nodes(:,2),1:3);     % second,
p3 = nodes2coord(faces2nodes(:,3),1:3);     % and third points of faces
vec1 = p2-p1;                               % two vectors defining
vec2 = p3-p1;                               % the face plane
normals = cross(vec1,vec2,2);               % normals

% outward normals of the element's four faces
% (note that there is no need to normalize the normals)
B_K_invT = amt(aminv(B_K));
n1 = squeeze(amsv(B_K_invT, [ 0  0 -1]))';
n2 = squeeze(amsv(B_K_invT, [ 0 -1  0]))';
n3 = squeeze(amsv(B_K_invT, [-1  0  0]))';
n4 = squeeze(amsv(B_K_invT, [ 1  1  1]))';

% calculate the product of the positive normals and the
% outward normals of the element faces
tmp(:,1) = sum(normals(elems2faces(:,1),:) .* n1 , 2);
tmp(:,2) = sum(normals(elems2faces(:,2),:) .* n2 , 2);
tmp(:,3) = sum(normals(elems2faces(:,3),:) .* n3 , 2);
tmp(:,4) = sum(normals(elems2faces(:,4),:) .* n4 , 2);

% if the outward normal is in the same direction as the
% chosen positive normal direction, it has a sign +1
signs = tmp ./ abs(tmp);

end