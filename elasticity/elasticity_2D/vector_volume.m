%
function f_V=vector_volume(ELEM,COORD,f_V_int,HatP,WEIGHT)
 
% ======================================================================
%
% Assembling of the vector of volume forces
%
%    output: 
%      f_V - vector of volume forces, size(f_V)=(2,n_n) where n_n is 
%            the number of nodes
%
%    input data:
%      ELEM    - to indicate nodes belonging to each element 
%                size(ELEM)=(n_p,n_e) where n_e is a number of elements 
%                and n_p is a number of the nodes within one element
%      COORD   - coordinates of the nodes, size(COORD)=(2,n_n)
%      f_V_int - values of volume forces at integration points
%                size(f_V_int)=(2,n_int), where n_int=n_e*n_q is a number
%                of integration points, n_q is a number of quadrature points
%      HatP    - values of the basis functions at quadrature points
%                size(HatP))=(n_p,n_q)
%      WEIGHT  - weight coefficients, size(WEIGHT)=(1,n_int)
%
% =========================================================================

%
% Auxilliary notation
%
  n_n=size(COORD,2);    % number of nodes
  n_e=size(ELEM,2);     % number of elements
  n_p=size(ELEM,1);     % number of vertices per element
  n_q=size(HatP,2);     % number of quadratic points
  n_int = n_e*n_q ;     % total number of integrations points
    
  % extension of the input array HatP by replication
  % size(HatPhi)=(n_p,n_int)
  HatPhi=repmat(HatP,1,n_e);

%
% Assembling of the vector of volume forces, size(f_V)=(2,n_n)
%
  % values at integration points, size(vF1)=size(vF2)=(n_p,n_int)   
  vF1 = HatPhi.*(ones(n_p,1)*(WEIGHT.*f_V_int(1,:)));    
  vF2 = HatPhi.*(ones(n_p,1)*(WEIGHT.*f_V_int(2,:)));   
  % row and column indices, size(iF)=size(jF)=(n_p,n_int)   
  iF = ones(n_p,n_int);
  jF = kron(ELEM,ones(1,n_q));
  % the asssembling by using the sparse command - values v for duplicate
  % doubles i,j are automatically added together
  f_V = [ sparse(iF(:), jF(:), vF1(:), 1, n_n);
          sparse(iF(:), jF(:), vF2(:), 1, n_n) ];
end      
  
