%
function f_V=vector_volume(ELEM,COORD,f_V_int,HatP,DHatP1,DHatP2,WF)
 
% ======================================================================
%
% Assembling of the vector of volume forces
%
%    output: f_V - vector of volume forces, size(f_V)=(2,n_n)
%                  where n_n is the number of nodes
%    input data:
%      ELEM - to indicate nodes belonging to each element 
%             size(ELEM)=(n_e,n_p) where n_e is a number of elements 
%             and n_p is a number of the nodes within one element                           
%      COORD - coordinates of the nodes, size(COORD)=(2,n_n)
%      f_V_int - values of volume forces at integration points
%                size(f_V_int)=(2,n_int), where n_int=n_e*n_q is a number
%                of integration points, n_q is a number of quadrature pnts.
%      HatP - values of the basis functions at quadrature points
%      DHatP1- xi_1-derivatives of the basis functions at quadrature points
%      DHatP2- xi_2-derivatives of the basis functions at quadrature points
%              size(HatP)=size(DHatP1)=size(DHatP2)=(n_p,n_q)
%      WF - weight factors at quadrature points, size(WF)=(1,n_q)
%
% =========================================================================

%
% Auxilliary notation
%
  n_n=size(COORD,2);    % number of nodes including midpoints
  n_e=size(ELEM,2);     % number of elements
  n_p=size(ELEM,1);     % number of vertices per element
  n_q=length(WF);       % number of quadratic points
  n_int = n_e*n_q ;     % total number of integrations points
    
%
% Jacobians and their determinants at volume integration point
%   

  % extension of the input arrays HatP,DHatP1,DHatP2 by replication
  % size(DHatPhi1)=size(DHatPhi2)=size(HatPhi)=(n_p, n_int)
  DHatPhi1=repmat(DHatP1,1,n_e);
  DHatPhi2=repmat(DHatP2,1,n_e);
  HatPhi  =repmat(HatP,1,n_e);
  
  % coordinates of nodes defining each element
  % size(COORDe1)=size(COORDe2)=(n_p, n_e)
  COORDe1=reshape(COORD(1,ELEM(:)),n_p,n_e);
  COORDe2=reshape(COORD(2,ELEM(:)),n_p,n_e);
  
  % coordinates of nodes around each integration point
  % size(COORDint1)=size(COORDint2)=(n_p, n_int)
  COORDint1=kron(COORDe1,ones(1,n_q)); 
  COORDint2=kron(COORDe2,ones(1,n_q));
  
  % components of the Jacobians: size=(1, n_int)
  J11=sum(COORDint1.*DHatPhi1); J12=sum(COORDint2.*DHatPhi1); 
  J21=sum(COORDint1.*DHatPhi2); J22=sum(COORDint2.*DHatPhi2); 
  
  % determinant of the Jacobian: size=(1, n_int)
  DET=abs(J11.*J22 - J12.*J21);
  
  % weight coefficients: size(WEIGHT)=(1, n_int)
  WEIGHT = DET.*repmat(WF,1,n_e);

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
  
