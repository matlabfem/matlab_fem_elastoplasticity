function [nodes2coord,elems2nodes,dirichlet]=refinement_uniform_2D(nodes2coord,elems2nodes,dirichlet) 
%function: [nodes2coord,elems2nodes,dirichlet]=refinement_uniform_2D(nodes2coord,elems2nodes,dirichlet) 
%requires: get_edges
%uniform refinement of a 2D triangulation 

%uniform refinement   
[elems2edges, edges2nodes]=get_edges(elems2nodes);    
nodes2edge=sparse(edges2nodes(:,1),edges2nodes(:,2),1:size(edges2nodes,1),size(nodes2coord,1),size(nodes2coord,1));
nodes2edge=symetrizeMatrix(nodes2edge); 
    
%elements on uniformly refined mesh
elements_internal=elems2edges+size(nodes2coord,1);
elements_refin1= [elems2nodes(:,1) elements_internal(:,3) elements_internal(:,2)];
elements_refin2= [elems2nodes(:,2) elements_internal(:,1) elements_internal(:,3)];
elements_refin3= [elems2nodes(:,3) elements_internal(:,2) elements_internal(:,1)];    
elems2nodes=[elements_internal; elements_refin1; elements_refin2; elements_refin3];  

if (nargin==3)
%dirichlet edges of uniformly refined mesh
dirichlet_edges=diag(nodes2edge(dirichlet(:,1),dirichlet(:,2)));
dirichlet=[dirichlet(:,1) dirichlet_edges+size(nodes2coord,1); dirichlet_edges+size(nodes2coord,1) dirichlet(:,2)];
end


%coordinates of uniformly refined mesh
coordinates_internal=(nodes2coord(edges2nodes(:,1),:)+nodes2coord(edges2nodes(:,2),:))/2;
nodes2coord=[nodes2coord; coordinates_internal];   

function A_sym = symetrizeMatrix(A)
[i,j,k]=find(A);
W=sparse([i; j], [j; i], ones(size(k,1)*2,1));
A_help=sparse([i; j], [j; i], [k; k]);
[i,j,k]=find(A_help);
[i,j,kk]=find(W);
A_sym=sparse(i,j,(kk.^(-1)).*k); %Now Kantennr_sym is a symetric form of Kantennr