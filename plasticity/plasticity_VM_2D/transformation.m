%
function Q_node=transformation(Q_int,elem,weight)
 
% ======================================================================
%
% Transformation of function values at integration points to function
% values at nodes of the finite element mesh.
%
%    output: 
%      Q_node - values of a function Q at nodes of the FE mesh,
%               size(Q_node)=(1,n_n), where n_n is the number of nodes
%                  
%    input data:
%      Q_int  - values of a function Q at integration points, size(Q_int)=(1,n_int) 
%      elem   - to indicate nodes belonging to each element 
%               size(elem)=(n_p,n_e) where n_e is a number of elements 
%               and n_p is a number of the nodes within one element                           
%      weight - weight factors at integration points, size(weight)=(1,n_int)
%
% =========================================================================

%
% Auxilliary notation
%
  n_e=size(elem,2);     % number of elements
  n_p=size(elem,1);     % number of vertices per element
  n_int=length(weight); % total number of integrations points
  n_q=n_int/n_e;        % number of quadrature points on one element 

%
% F1 - 1*n_n array, to each node we compute numerically the integral of Q 
%                   through a vicinity of the node
% F2 - 1*n_n array, to each node we compute the area of the vicinity
%

  % values at integration points, size(vF1)=size(vF2)=(n_p,n_int)   
  vF1 = ones(n_p,1)*(weight.*Q_int);    
  vF2 = ones(n_p,1)*weight; 

  % row and column indices, size(iF)=size(jF)=(n_p,n_int)   
  iF = ones(n_p,n_int);
  jF = kron(elem,ones(1,n_q));
  
  % the asssembling by using the sparse command - values v for duplicate
  % doubles i,j are automatically added together
  F1 = sparse(iF(:), jF(:), vF1(:)) ;
  F2 = sparse(iF(:), jF(:), vF2(:)) ;
  
%
% Approximated values of the function Q at nodes of the FE mesh
%
  Q_node = F1./F2;
  
