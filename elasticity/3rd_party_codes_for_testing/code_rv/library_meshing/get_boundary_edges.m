function bedges2edges = get_boundary_edges(edges2nodes,bedges2nodes)

% GETBOUNDARYEDGES
%   Calculates boundary edges for 2D and 3D meshes. This is needed for
%   setting boundary conditions for the eddy current problem.
%
% IN:  edges2nodes      edges by their nodes
%      bedges2nodes     boundary edges/faces by their nodes in 2D/3D
%
% OUT: bedges2edges     boundary edge indexes to 'edges2nodes'
%

dim = size(bedges2nodes,2);

if ( dim == 3 )
    %extracts sets of edges 
    edges1 = bedges2nodes(:,[1 2]);
    edges2 = bedges2nodes(:,[2 3]);
    edges3 = bedges2nodes(:,[3 1]);
    
    %as sets of their nodes (vertices)
    vertices = zeros(size(bedges2nodes,1)*3,2);
    vertices(1:3:end,:) = edges1;
    vertices(2:3:end,:) = edges2;
    vertices(3:3:end,:) = edges3;

    %repeated sets of nodes (joint edges) are eliminated 
    bedges2nodes = deleterepeatedrows(vertices);
end

matrix = [edges2nodes; bedges2nodes];
[matrixs,tags] = sortrows(sort(matrix,2));

% which ones were reps? k is a vector of indexes to matrixs.
k = find(all(diff(matrixs)==0,2));

% tags(k) is an index vector to edge2nodes (matrix) and denotes those edges
% which are on boundary
% tags(k+1) is an index vector to matrix and matrix(tags(k+a)) is the
% same as bedges2nodes, but in different order.

% we could just return tags(k), but we want that the order is the same
% as in bedges2nodes
[~,tags2]=sort(tags(k+1));
bedges2edges = tags(k(tags2));

end