function midp = get_midpoints( elems2nodes, nodes2coord )

%MIDPOINTS   
%   Calculate midpoints of the elements of the mesh
%
% SYNTAX:  midp = get_midpoints( nodes2coord )
%
% IN:   elems2nodes       elements by their nodes
%       nodes2coord       nodes by their coordinates
%       
%
% OUT:  midp              middle points of elements
%

dim = size(nodes2coord,2);

if ( dim==2 )
    midp = ( nodes2coord(elems2nodes(:,1),:) + ...
             nodes2coord(elems2nodes(:,2),:) + ...
             nodes2coord(elems2nodes(:,3),:) ) ./ 3;
elseif (dim==3)
    midp = ( nodes2coord(elems2nodes(:,1),:) + ...
             nodes2coord(elems2nodes(:,2),:) + ...
             nodes2coord(elems2nodes(:,3),:) + ...
             nodes2coord(elems2nodes(:,4),:) ) ./ 4;
else
    error('GET_MIDPOINTS: An error calculating the midpoints of elements.')
end

end

