function [nodes2coord,elems2nodes,dirichlet]=refinement_uniform_3D(nodes2coord,elems2nodes,dirichlet) 
%function: [nodes2coord,elems2nodes,dirichlet]=refinement_uniform_3D(nodes2coord,elems2nodes,dirichlet) 
%uniform refinement of a 3D triangulation 

%uniform refinement   
[nodes2coord,elems2nodes]=tetrarefine3(nodes2coord,elems2nodes);
%Dirichlet not working yet!!!

%[coordinates,elements3]=tetrarefine4(coordinates,elements3,[],[],[]);