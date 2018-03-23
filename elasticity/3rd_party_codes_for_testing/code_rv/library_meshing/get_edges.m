function [elems2edges, edges2nodes]=get_edges(elems2nodes)
%function: [element2edges, edge2nodes]=get_edges(elems2nodes)
%requires: deleterepeatedrows
%generates edges of (triangular) triangulation defined in elems2nodes
%elems2nodes is matrix, whose rows contain numbers of its element nodes 
%element2edges returns edges numbers of each triangular element
%edge2nodes returns two node numbers of each edge
%example in 2D: [element2edges, edge2nodes]=get_edges([1 2 3; 2 4 3])
%example in 3D: [element2edges, edge2nodes]=get_edges([1 2 3 4; 1 2 3 5; 1 2 4 6])

%2D case
if (size(elems2nodes,2)==3)
    %extracts sets of edges 
    edges1=elems2nodes(:,[2 3]);
    edges2=elems2nodes(:,[3 1]);
    edges3=elems2nodes(:,[1 2]);

    %as sets of their nodes (vertices)
    vertices=zeros(size(elems2nodes,1)*3,2);
    vertices(1:3:end,:)=edges1;
    vertices(2:3:end,:)=edges2;
    vertices(3:3:end,:)=edges3;

    %repeated sets of nodes (joint edges) are eliminated 
    [edges2nodes,elems2edges]=deleterepeatedrows(vertices);
    elems2edges=reshape(elems2edges,3,size(elems2nodes,1))';
end

%3D case
if (size(elems2nodes,2)==4)
    %extracts sets of edges 
    edges1=elems2nodes(:,[1 2]);
    edges2=elems2nodes(:,[1 3]);
    edges3=elems2nodes(:,[1 4]);
    edges4=elems2nodes(:,[2 3]);
    edges5=elems2nodes(:,[3 4]);
    edges6=elems2nodes(:,[4 2]);
    
    %as sets of their nodes (vertices)
    vertices=zeros(size(elems2nodes,1)*6,2);
    vertices(1:6:end,:)=edges1;
    vertices(2:6:end,:)=edges2;
    vertices(3:6:end,:)=edges3;
    vertices(4:6:end,:)=edges4;
    vertices(5:6:end,:)=edges5;
    vertices(6:6:end,:)=edges6;
    
    %repeated sets of nodes (joint edges) are eliminated 
    [edges2nodes,elems2edges]=deleterepeatedrows(vertices);
    elems2edges=reshape(elems2edges,6,size(elems2nodes,1))';
end
