function [coordinates,elements,dirichlet,neumann]=refinement_uniform(coordinates,elements,dirichlet,neumann) 
%function: [coordinates,elements3,dirichlet]=refinement_uniform(coordinates,elements3,dirichlet) 
%requires: getEdges
%uniform refinement of a 2D triangulation 

%uniform refinement   
[edge2nodes, ~, element2edges, ~]=getEdges(elements);
nodes2edge_sparse=sparse(edge2nodes(:,1),edge2nodes(:,2),1:size(edge2nodes,1),size(coordinates,1),size(coordinates,1));
nodes2edge_sparse=symetrizeMatrix(nodes2edge_sparse); 
    
%elements on uniformly refined mesh
elements3_internal=element2edges+size(coordinates,1);
elements3_refin1= [elements(:,1) elements3_internal(:,3) elements3_internal(:,2)];
elements3_refin2= [elements(:,2) elements3_internal(:,1) elements3_internal(:,3)];
elements3_refin3= [elements(:,3) elements3_internal(:,2) elements3_internal(:,1)];    
elements=[elements3_internal; elements3_refin1; elements3_refin2; elements3_refin3];  

if (nargin>=3)
    %dirichlet edges of uniformly refined mesh
    dirichlet_edges=diag(nodes2edge_sparse(dirichlet(:,1),dirichlet(:,2)));
    dirichlet=[dirichlet(:,1) dirichlet_edges+size(coordinates,1); dirichlet_edges+size(coordinates,1) dirichlet(:,2)];
end

if (nargin==4)
    if ~isempty(neumann)
        %neumann edges of uniformly refined mesh
        neumann_edges=diag(nodes2edge_sparse(neumann(:,1),neumann(:,2)));
        neumann=[neumann(:,1) neumann_edges+size(coordinates,1); neumann_edges+size(coordinates,1) neumann(:,2)];
    end
end

%coordinates of uniformly refined mesh
coordinates_internal=(coordinates(edge2nodes(:,1),:)+coordinates(edge2nodes(:,2),:))/2;
coordinates=[coordinates; coordinates_internal];   

function A_sym = symetrizeMatrix(A)
[i,j,k]=find(A);
W=sparse([i; j], [j; i], ones(size(k,1)*2,1));
A_help=sparse([i; j], [j; i], [k; k]);
[i,j,k]=find(A_help);
[i,j,kk]=find(W);
A_sym=sparse(i,j,(kk.^(-1)).*k); %Now Kantennr_sym is a symetric form of Kantennr