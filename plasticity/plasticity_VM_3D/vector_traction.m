function f_t = vector_traction(ELEM_s,COORD,f_t_int,HatP_s,DHatP1_s,DHatP2_s,WF_s)

% =========================================================================
%
% Assembling of the vector of traction forces acting on the upper side of
% the 3D body
%
%    output: 
%      f_t - vector of traction forces, size(f_V)=(2,n_n), where n_n is 
%            the number of nodes
%
%    input data:
%      ELEM_s   - to indicate nodes belonging to each surface element 
%                 size(ELEM_s)=(n_p_s,n_e_s) where n_e_s is a number of 
%                 surface elements n_p_s is a number of the nodes within one
%                 surface element                           
%      COORD    - coordinates of the nodes, size(COORD)=(3,n_n)
%      f_t_int  - values of traction forces at integration points
%                 size(f_t_int)=(3,n_int_s), where n_int_s=n_e_s*n_q_s is a 
%                 number of surface integration points and n_q_s is a number  
%                 of quadrature points on the surface
%      HatP_s   - values of the surface basis functions at quadrature points
%      DHatP1_s - xi_1-derivatives of the surface basis functions at q. p.
%      DHatP2_s - xi_2-derivatives of the surface basis functions at q. p.
%                 size(HatP_s)=size(DHatP1_s)=size(DHatP2_s)=(n_p_s,n_q_s)
%      WF_s     - weight factors at surface quadrature points, 
%                 size(WF)=(1,n_q_s)
%
% =========================================================================

%
% Auxilliary notation
%

  n_n=size(COORD,2);    % number of nodes including midpoints
  n_e_s=size(ELEM_s,2); % number of surface elements
  n_p_s=size(ELEM_s,1); % number of nodes within one surface element
  n_q_s=length(WF_s);   % number of quadrature points on a surface element
  n_int_s=n_e_s*n_q_s ; % number of integration points on the surface
                        % (on the upper side of the body)

%
% Jacobians and their determinants at surface integration points
%                         

  % extension of the input arrays HatP_s,DHatP1_s,DHatP2_s,DHatP3_s by
  % replication, size(DHatPhi1_s)=size(DHatPhi2_s)=(n_p_s,n_int_s)
  DHatPhi1_s=repmat(DHatP1_s,1,n_e_s);
  DHatPhi2_s=repmat(DHatP2_s,1,n_e_s);
    HatPhi_s=repmat(HatP_s,1,n_e_s)  ;
  
  % coordinates of nodes defining each surface element
  % size(COORDs1)=size(COORDs2)=size(COORDs3)=(n_p_s,n_e_s)
  COORDs1=reshape(COORD(1,ELEM_s(:)),n_p_s,n_e_s);
  COORDs2=reshape(COORD(3,ELEM_s(:)),n_p_s,n_e_s);
  
  % coordinates of nodes around each surface integration point
  % size(COORDint1)=size(COORDint2)=(n_p_s n_int_s)
  COORDint1=kron(COORDs1,ones(1,n_q_s)); 
  COORDint2=kron(COORDs2,ones(1,n_q_s));
  
  % components of the Jacobians: size=(1,n_int_s)
  J11=sum(COORDint1.*DHatPhi1_s); J12=sum(COORDint2.*DHatPhi1_s); 
  J21=sum(COORDint1.*DHatPhi2_s); J22=sum(COORDint2.*DHatPhi2_s);
  
  % determinant of the Jacobian: size=(1,n_int_s)
  DET_s=J11.*J22-J12.*J21;
  
  % weight coefficients: size(WEIGHT)=(1,n_int_s)
  WEIGHT_s = abs(DET_s).*repmat(WF_s,1,n_e_s);
  
%
% Assembling of the vector of traction forces, size(f_V)=(3,n_n)
%

  % auxilliary values at surface integration points, 
  % size(vF1)=size(vF2)=size(vF2)=(n_p_s,n_int_s)   
  vF1 = HatPhi_s.*(ones(n_p_s,1)*(WEIGHT_s.*f_t_int(1,:)));    
  vF2 = HatPhi_s.*(ones(n_p_s,1)*(WEIGHT_s.*f_t_int(2,:)));   
  vF3 = HatPhi_s.*(ones(n_p_s,1)*(WEIGHT_s.*f_t_int(3,:)));   
  % row and column indices, size(iF)=size(jF)=(n_p_s,n_int_s)   
  iF = ones(n_p_s,n_int_s);
  jF = kron(ELEM_s,ones(1,n_q_s));
  % the asssembling by using the sparse command - values v for duplicate
  % doubles i,j are automatically added together
  f_t = [ sparse(iF(:), jF(:), vF1(:), 1, n_n);
          sparse(iF(:), jF(:), vF2(:), 1, n_n);
          sparse(iF(:), jF(:), vF3(:), 1, n_n) ];  
  
 end