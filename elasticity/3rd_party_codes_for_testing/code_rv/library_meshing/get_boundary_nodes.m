function bnodes = get_boundary_nodes(elems2nodes)
%works only on 2D triangular and 3D tetrahedral meshes!

%2D case
if (size(elems2nodes,2)==3)
    if 0
        [elems2edges, edges2nodes]=get_edges(elems2nodes)   
        %add on
        edge2elems=entryInWhichRows(elems2edges); 
        if (size(edge2elems,2)==1)
            edge2elems=[edge2elems 0*edge2elems];    %all edges are boundary edges!!!
        end 
        node2edges=entryInWhichRows(edge2nodes); 
    else
        [edge2nodes, edge2elems, ~, ~]=getEdges(elems2nodes);
    end
    
    I= edge2elems(:,2)==0;
    bnodes=unique(edge2nodes(I,:));
end

%3D case
if (size(elems2nodes,2)==4)
    [elems2faces, faces2nodes]=get_faces(elems2nodes);
    
    %add on
    face2elems=entryInWhichRows(elems2faces); 
    if (size(face2elems,2)==1)
        face2elems=[face2elems 0*face2elems];    %all edges are boundary edges!!!
    end 
    node2faces=entryInWhichRows(faces2nodes); 
    
    I= face2elems(:,2)==0;
    bnodes=unique(faces2nodes(I,:));
    
end
