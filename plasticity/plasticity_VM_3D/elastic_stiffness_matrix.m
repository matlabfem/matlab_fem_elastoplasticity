
function  [K,B,WEIGHT,iD,jD,D]=elastic_stiffness_matrix(ELEM,COORD,...
                                       shear,bulk,DHatP1,DHatP2,DHatP3,WF)
% =========================================================================
%
%  This function assembles the elastic stiffness matrix and some other
%  auxilliary arrays which are important for assembling of the tangent
%  stiffness matrix related to the elastoplacitic body
%
%  input data:
%    ELEM   - array containing numbers of nodes defining each element,
%             size(ELEM)=(n_p,n_e), n_e = number of elements
%    COORD  - coordinates of the nodes, size(coord)=(3,n_n) where n_n is a
%             number of nodes
%    shear  - shear moduli at integration points, size(shear)=(1,n_int)
%    bulk   - bulk moduli at integration points, size(bulk)=(1,n_int)
%    DHatP1 - derivatives of basis functions at the quadrature points 
%             in the direction xi_1, size(HatP)=(n_p,n_q)
%    DHatP2 - derivatives of basis functions at the quadrature points 
%             in the direction xi_2, size(DHatP2)=(n_p,n_q)
%    DHatP3 - derivatives of basis functions at the quadrature points 
%             in the direction xi_3, size(DHatP3)=(n_p,n_q) 
%    WF     - weight factors on the reference element, size(WF)=(1,n_q)
%
%  output data:
%    K      - the elastic stiffness matrix, size(K)=(3*n_n, 3*n_n)
%    B      - the strain-displacement matrix, size(B)=(6*n_int,3*n_n)
%    WEIGHT - the weight coefficients for each quadrature point, 
%             size(WEIGHT)=(1,n_int)
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
  
%
% Jacobian, its determinant and inverse, 
% derivative of local basis functions
% 

  % extension of the input arrays DHatP1,DHatP2,DHatP3 by replication
  % size(DHatPhi1)=size(DHatPhi2)=size(DHatPhi3)=(n_p,n_int)
  DHatPhi1=repmat(DHatP1,1,n_e);
  DHatPhi2=repmat(DHatP2,1,n_e);
  DHatPhi3=repmat(DHatP3,1,n_e);
  
  % coordinates of nodes defining each element
  % size(COORDe1)=size(COORDe2)=size(COORDe3)=(n_p,n_e)
  COORDe1=reshape(COORD(1,ELEM(:)),n_p,n_e);
  COORDe2=reshape(COORD(2,ELEM(:)),n_p,n_e);
  COORDe3=reshape(COORD(3,ELEM(:)),n_p,n_e);
  
  % coordinates of nodes around each integration point
  % size(COORDint1)=size(COORDint2)=size(COORDint3)=(n_p,n_int)
  COORDint1=kron(COORDe1,ones(1,n_q)); 
  COORDint2=kron(COORDe2,ones(1,n_q));
  COORDint3=kron(COORDe3,ones(1,n_q));  
  
  % components of the Jacobians: size=(1,n_int)
  J11=sum(COORDint1.*DHatPhi1); J12=sum(COORDint2.*DHatPhi1); J13=sum(COORDint3.*DHatPhi1);
  J21=sum(COORDint1.*DHatPhi2); J22=sum(COORDint2.*DHatPhi2); J23=sum(COORDint3.*DHatPhi2);
  J31=sum(COORDint1.*DHatPhi3); J32=sum(COORDint2.*DHatPhi3); J33=sum(COORDint3.*DHatPhi3);
  
  % determinant of the Jacobian: size=(1,n_int)
  DET=J11.*(J22.*J33-J32.*J23) - J12.*(J21.*J33-J23.*J31) + J13.*(J21.*J32-J22.*J31);
  
  % components of the inverse to the Jacobian: size=(1,n_int)
  Jinv11 =  (J22.*J33-J23.*J32)./DET; Jinv12 = -(J12.*J33-J13.*J32)./DET; Jinv13 =  (J12.*J23-J13.*J22)./DET; 
  Jinv21 = -(J21.*J33-J23.*J31)./DET; Jinv22 =  (J11.*J33-J13.*J31)./DET; Jinv23 = -(J11.*J23-J13.*J21)./DET; 
  Jinv31 =  (J21.*J32-J22.*J31)./DET; Jinv32 = -(J11.*J32-J12.*J31)./DET; Jinv33 =  (J11.*J22-J12.*J21)./DET; 
  
  % derivatives of local basis functions w.r.t the coordinates x_1,x_2,x_3:
  % size(DPhi1)=size(DPhi2)=size(DPhi3)=(n_p,n_int)
  DPhi1 = repmat(Jinv11,n_p,1).*DHatPhi1 + repmat(Jinv12,n_p,1).*DHatPhi2 + repmat(Jinv13,n_p,1).*DHatPhi3;
  DPhi2 = repmat(Jinv21,n_p,1).*DHatPhi1 + repmat(Jinv22,n_p,1).*DHatPhi2 + repmat(Jinv23,n_p,1).*DHatPhi3;
  DPhi3 = repmat(Jinv31,n_p,1).*DHatPhi1 + repmat(Jinv32,n_p,1).*DHatPhi2 + repmat(Jinv33,n_p,1).*DHatPhi3;
  
%
% assembling of the strain-displacement matrix
% size(B)=(6*n_int,3*n_n)
%   
  
  % values of the strain-displacement matrix B
  n_b = 18*n_p ;
  vB=zeros(n_b,n_int);        % size(vB)=(18*n_p,n_int)
  vB(1:18:n_b-17,:)=DPhi1; vB(10:18:n_b- 8,:)=DPhi1; vB(18:18:n_b  ,:)=DPhi1;
  vB(4:18:n_b-14,:)=DPhi2; vB( 8:18:n_b-10,:)=DPhi2; vB(17:18:n_b-1,:)=DPhi2;
  vB(6:18:n_b-12,:)=DPhi3; vB(11:18:n_b- 7,:)=DPhi3; vB(15:18:n_b-3,:)=DPhi3;

  % i-th and j-th indices of B: size(iB)=size(jB)=(18*n_p,n_int)
  AUX=reshape(1:6*n_int,6,n_int);
  iB=repmat(AUX,3*n_p,1);  
  
  AUX1 = [1;1;1]*(1:n_p);     % size(AUX1)=(3,n_p)
  AUX2 = [2;1;0]*ones(1,n_p); % size(AUX2)=(3,n_p)
  AUX3=3*ELEM((AUX1(:))',:)-kron(ones(1,n_e),AUX2(:));
                              % size(AUX3)=(3*n_p,n_p)
  jB=kron(AUX3,ones(6,n_q));
  
  % the sparse strain-displacement matrix B
  B=sparse(iB(:),jB(:),vB(:), 6*n_int,3*n_n);  

%  
% assembling of the elastic stress-strain matrix 
% size(D)=(6*n_int,6*n_int)
%
  
  % elastic tensor at integration points:
  IOTA=[1;1;1;0;0;0];  
  VOL=IOTA*IOTA'; 
  DEV=diag([1,1,1,1/2,1/2,1/2])-VOL/3; 
  ELAST=2*DEV(:)*shear+VOL(:)*bulk;   % size(ELAST)=(36,n_int)
    
  % weight coefficients: size(WEIGHT)=(1,n_int)
  WEIGHT = abs(DET).*repmat(WF,1,n_e);

  % assembling of the sparse matrix D
  iD=repmat(AUX,6,1); 
  jD=kron(AUX,ones(6,1));
  vD=ELAST.*(ones(36,1)*WEIGHT);
  D=sparse(iD,jD,vD);

%
% elastic stiffness matrix: size(K)=(3*n_n,3*n_n)
%
  K = B'*D*B ;    
 
end  % end of function