
function [S,IND_p,DS_p,Ep,Hard]=constitutive_problem(E,Ep_prev,Hard_prev,...
                                                      shear,bulk,a,Y)

% =========================================================================
%
% The aim of this function is to construct constitutive and consistent 
% tangent operators at integration points 1,2,...,n_int. These operators
% are related to elastoplastic model containing the von Mises yield
% criterion and a linear kinematic hardening.
%
% Input data:
%  E_new - current strain tensor, size(E_new)=(6, n_int)
%  Ep_prev - plastic strain tensor from the previous time step,
%            size(Ep_prev)=(6, n_int)
%  Hard_prev - kinematic hardening (tensor) from the previous time step
%              size(Hard_prev)=(6, n_int)
%  ELAST - elastic fourth order tensor at integration points
%          size(ELAST)=(36, n_int)
%  DEV - deviatoric operator, size(DEV)=(6, 6)
%  shear - shear moduli at integration points, size(shear)=(1, n_int)
%  bulk - bulk moduli at integration points, size(shear)=(1, n_int)
%  a - hardening parameters at integration points, size(a)=(1, n_int)
%  Y - initial yield stress at integration points, size(Y)=(1, n_int)
%
% Output data:
%  S - stress tensors at integration points, size(S)=(6, n_int)
%  IND_p - 1*n_int logical array indicating integration points with plastic 
%          response, n_plast=number of the points with plastic response
%  DS_p - plastic parts of consistent tangent operators at integr. points,
%         size(D_p)=(36, n_plast)
%  Ep - plastic strain, size(Ep)=(6,n_int)
%  Hard - tensors of kinematic hardening, size(Hard)=(6, n_int)
%
% =========================================================================
%

%
% Elastic tensor at integration points, size(ELAST)=(36, n_int).
% Deviatoric and volumetric 6x6 matrices
%
  IOTA=[1;1;1;0;0;0];  
  VOL=IOTA*IOTA'; 
  DEV=diag([1,1,1,1/2,1/2,1/2])-VOL/3; 
  ELAST=2*DEV(:)*shear+VOL(:)*bulk; 

%
% Trial variables:
%   E_tr - trial strain tensors, size(E_tr)=(6,n_int)
%   S_tr - trial stress tensors, size(E_tr)=(6,n_int)
%   SD_tr - deviatoric part of the diffence between trial strain 
%          and hardening tensors, size(SD_tr)=(6,n_int)
%   norm_SD - norm of SD_tr, size(norm_SD)=(1,n_int)
%

  E_tr=E-Ep_prev;                          
  S_tr=2*repmat(shear,6,1).*(DEV*E_tr)+repmat(bulk,6,1).*(VOL*E_tr);   
  SD_tr=DEV*(2*repmat(shear,6,1).*E_tr)-Hard_prev;  
  norm_SD=sqrt(sum(SD_tr(1:3,:).*SD_tr(1:3,:))+2*sum(SD_tr(4:6,:).*SD_tr(4:6,:)));  
  
%
% Evaluation of the yield criterion and specification of integration points
% with plastic response
%

  CRIT=norm_SD-Y;  % size(CRIT)=(1, n_int)
  IND_p=CRIT>0; % logical array, size(INDplast)=(1, n_int)   

%
% The elastic prediction of unknowns
%

  S=S_tr; 

%
% The plastic correction at the selected integration points
%

  % N_hat - normalized tensor to SD_tr, size(N_hat)=(6, n_int)
  N_hat=SD_tr(:,IND_p)./repmat(norm_SD(IND_p),6,1);
  
  % plastic multipliers at integration points with plastic response 
  denom = 2*shear(IND_p)+a(IND_p);  
  lambda=CRIT(IND_p)./denom; 
  
  % correction of the stress tensors
  S(:,IND_p)=S(:,IND_p)-repmat(2*shear(IND_p).*lambda,6,1).*N_hat;
 
  % correction of the consistent tangent operator
  ID=DEV(:)*ones(1,length(denom));
  NN_hat=repmat(N_hat,6,1).*kron(N_hat,ones(6,1));
  const=((2*shear(IND_p)).^2)./denom;
  DS_p=-repmat(const,36,1).*ID+...
         repmat((const.*Y(IND_p))./norm_SD(IND_p),36,1).*(ID-NN_hat);

%     
% Update of the plastic strain and the hardening variable (optional)
%

  if nargout>3    
     % plastic strain
     Ep=Ep_prev;
     Ep(:,IND_p)=Ep(:,IND_p)+([1;1;1;2;2;2]*lambda).*N_hat;
     % hardening 
     Hard=Hard_prev;
     Hard(:,IND_p)=Hard(:,IND_p)+repmat(a(IND_p).*lambda,6,1).*N_hat;
  end

 end
