
function  [K,B,WEIGHT,iD,jD,D]=elastic_stiffness_matrix(ELEM,COORD,...
                                            shear,bulk,DHatP1,DHatP2,WF)
 
% =========================================================================
%
%  This function assembles the elastic stiffness matrix and some other
%  auxilliary arrays which are important for assembling of the tangent
%  stiffness matrix related to the elastoplacitic body
%
%  input data:
%    ELEM - n_p x n_e array containing numbers of nodes defining each
%           element, n_e = number of elements
%    COORD - coordinates of the nodes, size(coord)=(2,n_n) where n_n is a
%            number of nodes
%    shear, bulk - material parameters at each integration point, 
%                  size(shear)=size(bulk)=(1,n_int)
%    DHatP1 - derivatives of basis functions at the quadrature points 
%             in the direction xi_1, size(HatP)=(n_p,n_q)
%    DHatP2 - derivatives of basis functions at the quadrature points 
%             in the direction xi_2, size(DHatP2)=(n_p,n_q)
%        WF - weight factors on the reference element, size(WF)=(1,n_q)
%
%  output data:
%    K - the elastic stiffness matrix, size(K)=(2*n_n, 2*n_n)
%    B - the strain-displacement matrix, size(B)=(3*n_int,2*n_n)
%    WEIGHT - the weight coefficients for each quadrature point, 
%             size(WEIGHT)=(1, n_int)
%
% ======================================================================
%

%
% auxilliary notation
%

  n_n=size(COORD,2);    % number of nodes including midpoints
  n_e=size(ELEM,2);     % number of elements
  n_p=size(ELEM,1);     % number of vertices per element
  n_q=length(WF);       % number of quadrature points
  n_int = n_e*n_q ;     % total number of integrations points
  
  if numel(shear)==1    %homogeneous material
     shear=shear*ones(1,n_int);
  end
  if numel(bulk)==1
     bulk=bulk*ones(1,n_int);
  end
     
%
% Jacobian, its determinant and inverse, 
% derivative of local basis functions
% 

  % extension of the input arrays DHatP1,DHatP2 by replication
  % size(DHatPhi1)=size(DHatPhi2)=(n_p, n_int)
  DHatPhi1=repmat(DHatP1,1,n_e);
  DHatPhi2=repmat(DHatP2,1,n_e);
  
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
  DET=J11.*J22 - J12.*J21;
  
  % components of the inverse to the Jacobian: size=(1, n_int)
  Jinv11 =  J22./DET; Jinv12 = -J12./DET;
  Jinv21 = -J21./DET; Jinv22 =  J11./DET; 
  
  DET=abs(DET);
  
  % derivatives of local basis functions w.r.t the coordinates x_1,x_2:
  % size(DPhi1)=size(DPhi2)=(n_p, n_int)
  DPhi1 = repmat(Jinv11,n_p,1).*DHatPhi1 + repmat(Jinv12,n_p,1).*DHatPhi2;
  DPhi2 = repmat(Jinv21,n_p,1).*DHatPhi1 + repmat(Jinv22,n_p,1).*DHatPhi2;
  
%
% assembling of the strain-displacement matrix
% size(B)=(3*n_int, 2*n_n)
%   
  
  % values of the strain-displacement matrix B
  n_b = 6*n_p ;
  vB=zeros(n_b,n_int);        % size(vB)=(6*n_p,n_int)
  vB(1:6:n_b-5,:)=DPhi1; vB(6:6:n_b  ,:)=DPhi1; 
  vB(5:6:n_b-1,:)=DPhi2; vB(3:6:n_b-3,:)=DPhi2; 

  % i-th and j-th indices of B: size(iB)=size(jB)=(6*n_p,n_int)
  AUX=reshape(1:3*n_int,3,n_int);
  iB=repmat(AUX,2*n_p,1);  
  
  AUX1 = [1;1]*(1:n_p);     % size(AUX1)=(2,n_p)
  AUX2 = [1;0]*ones(1,n_p); % size(AUX2)=(2,n_p)
  AUX3=2*ELEM((AUX1(:))',:)-kron(ones(1,n_e),AUX2(:));
                              % size(AUX3)=(2*n_p,n_p)
  jB=kron(AUX3,ones(3,n_q));
  
  % the sparse strain-displacement matrix B
  B=sparse(iB(:),jB(:),vB(:), 3*n_int,2*n_n);  

%  
% assembling of the elastic stress-strain matrix 
% size(D)=(3*n_int, 3*n_int)
%
  
  % elastic tensor at each integration point: 
  IOTA=[1;1;0];  
  VOL=IOTA*IOTA'; 
  IDENT=diag([1,1,1/2]); 
  DEV=IDENT-VOL/3; 
  ELAST=2*DEV(:)*shear+VOL(:)*bulk; % size(ELAST)=(9, n_int)
 
  % weight coefficients: size(WEIGHT)=(1, n_int)
  WEIGHT = DET.*repmat(WF,1,n_e);
  
  % assemblinng of the sparse matrix D
  iD=repmat(AUX,3,1); 
  jD=kron(AUX,ones(3,1));
  vD=ELAST.*(ones(9,1)*WEIGHT);
  D=sparse(iD,jD,vD);

%
% elastic stiffness matrix: size(K)=(2*n_n, 2*n_n)
%
  K = B'*D*B ;    
 
end  % end of function