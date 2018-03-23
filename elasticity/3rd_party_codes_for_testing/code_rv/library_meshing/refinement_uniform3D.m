function [coordinates,elements3,dirichlet]=refinement_uniform3D(coordinates,elements3,dirichlet) 
%function: [coordinates,elements3,dirichlet]=refinement_uniform(coordinates,elements3,dirichlet) 
%uniform refinement of a 3D triangulation 

%uniform refinement   
[coordinates,elements3]=tetrarefine3(coordinates,elements3);
%Dirichlet not working yet!!!

%[coordinates,elements3]=tetrarefine4(coordinates,elements3,[],[],[]);