
function [S,DS,IND_p,Ep,Hard]=constitutive_problem(E,Ep_prev,Hard_prev,...
                                                      shear,bulk,a,Y)

% =========================================================================
%
% The aim of this function is to construct constitutive and consistent 
% tangent operators at integration points 1,2,...,n_int. These operators
% are related to elastoplastic model containing the von Mises yield
% criterion and a linear kinematic hardening.
%
% Input data:
%  E         - current strain tensor, size(E)=(3,n_int)
%  Ep_prev   - plastic strain tensor from the previous time step,
%              size(Ep_prev)=(4,n_int)
%  Hard_prev - kinematic hardening (tensor) from the previous time step
%              size(Hard_prev)=(4,n_int)
%  shear     - shear moduli at integration points, size(shear)=(1,n_int)
%  bulk      - bulk moduli at integration points, size(shear)=(1,n_int)
%  a         - hardening parameters at integration points, size(a)=(1,n_int)
%  Y         - initial yield stress at integration points, size(Y)=(1,n_int)
%
% Output data:
%  S     - stress tensors at integration points, size(S)=(4, n_int)
%  DS    - consistent tangent operators at integr. points,
%          size(DS)=(9,n_plast)
%  IND_p - logical array indicating integration points with plastic response, 
%          size(IND_p)=(1,n_int), 
%          n_plast=number of the points with plastic response
%  Ep    - plastic strain, size(Ep)=(4,n_int)
%  Hard  - tensors of kinematic hardening, size(Hard)=(4,n_int)
%
% =========================================================================
%
  n_int=length(shear); % number of integration points

%
% Elastic tensor at integration points, size(ELAST)=(9,n_int).
% Deviatoric and volumetric 4x4 matrices
%
  IOTA=[1;1;0;1];  
  VOL=IOTA*IOTA'; 
  DEV=diag([1,1,1/2,1])-VOL/3; 
  Dev=DEV(1:3,1:3);
  Vol=VOL(1:3,1:3);

%
% Trial variables:
%   E_tr    - trial strain tensors, size(E_tr)=(4,n_int)
%   S_tr    - trial stress tensors, size(S_tr)=(4,n_int)
%   SD_tr   - deviatoric part of the diffence between trial strain 
%             and hardening tensors, size(SD_tr)=(4,n_int)
%   norm_SD - norm of SD_tr, size(norm_SD)=(1,n_int)
%

  E4=[E;zeros(1,n_int)];
  E_tr=E4-Ep_prev;                          
  S_tr=2*repmat(shear,4,1).*(DEV*E_tr)+repmat(bulk,4,1).*(VOL*E_tr);   
  SD_tr=DEV*(2*repmat(shear,4,1).*E_tr)-Hard_prev;  
  norm_SD=sqrt(sum(SD_tr([1,2,4],:).*SD_tr([1,2,4],:))+2*(SD_tr(3,:).*SD_tr(3,:)));  
  
%
% Evaluation of the yield criterion and specification of integration points
% with plastic response
%

  CRIT=norm_SD-Y;  % size(CRIT)=(1,n_int)
  IND_p=CRIT>0;    % logical array, size(IND_p)=(1,n_int)   

%
% The elastic prediction of unknowns
%

  S=S_tr; DS=2*Dev(:)*shear+Vol(:)*bulk;

%
% The plastic correction at the selected integration points
%

  % N_hat - normalized tensor to SD_tr, size(N_hat)=(4,n_int)
  N_hat=SD_tr(:,IND_p)./repmat(norm_SD(IND_p),4,1);
  
  % plastic multipliers at integration points with plastic response 
  denom = 2*shear(IND_p)+a(IND_p);  
  lambda=CRIT(IND_p)./denom; 
  
  % correction of the stress tensors
  S(:,IND_p)=S(:,IND_p)-repmat(2*shear(IND_p).*lambda,4,1).*N_hat;
 
  % correction of the consistent tangent operator
  ID=Dev(:)*ones(1,length(denom));
  NN_hat=repmat(N_hat(1:3,:),3,1).*kron(N_hat(1:3,:),ones(3,1));
  const=((2*shear(IND_p)).^2)./denom;
  DS(:,IND_p)=DS(:,IND_p)-repmat(const,9,1).*ID+...
         repmat((const.*Y(IND_p))./norm_SD(IND_p),9,1).*(ID-NN_hat);

%     
% Update of the plastic strain and the hardening variable (optional)
%

  if nargout>3    
     % plastic strain
     Ep=Ep_prev;
     Ep(:,IND_p)=Ep(:,IND_p)+([1;1;2;1]*lambda).*N_hat;
     % hardening 
     Hard=Hard_prev;
     Hard(:,IND_p)=Hard(:,IND_p)+repmat(a(IND_p).*lambda,4,1).*N_hat;
  end

 end
