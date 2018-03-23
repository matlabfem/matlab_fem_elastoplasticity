
function [S,DS,IND_p,lambda,Ep]=constitutive_problem(E,Ep_prev,shear,bulk,eta,c)

% =========================================================================
%
% The aim of this function is to construct constitutive and consistent 
% tangent operators at integration points 1,2,...,n_int. These operators
% are related to elastic-perfectly plastic model containing the 
% Drucker-Prager yield criterion.
%
% Input data:
%  E_new - current strain tensor, size(E_new)=(6, n_int)
%  Ep_prev - plastic strain tensor from the previous time step,
%            size(Ep_prev)=(6, n_int)
%  ELAST - elastic fourth order tensor at integration points
%          size(ELAST)=(36, n_int)
%  DEV - deviatoric operator, size(DEV)=(6, 6)
%  shear - shear moduli at integration points, size(shear)=(1, n_int)
%  bulk - bulk moduli at integration points, size(shear)=(1, n_int)
%  eta,c - inelastic parameters at integration points, 
%          size(eta)=size(c)=(1, n_int)
%
% Output data:
%  S - stress tensors at integration points, size(S)=(6, n_int)
%  IND_p - 1*n_int logical array indicating integration points with plastic 
%          response, n_plast=number of the points with plastic response
%  DS - consistent tangent operators at integr. points,
%         size(DS_p)=(36, n_plast)
%  lambda - plastic multipliers, size(lambda)=(1,n_int)
%  Ep - plastic strains, size(Ep)=(6,n_int)
%
% =========================================================================
%
  n_int=length(shear); % number of integration points

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
%   dev_E - deviatoric part of the trial strain, size(SD_tr)=(6,n_int)
%   norm_E - norm of dev_E, size(norm_SD)=(1,n_int)
%

  E_tr=E-Ep_prev;                          
  S_tr=2*repmat(shear,6,1).*(DEV*E_tr)+repmat(bulk,6,1).*(VOL*E_tr);   
  dev_E=DEV*E_tr;                       % deviatoric part of E_tr
  norm_E=sqrt(max(0,sum(E_tr.*dev_E))); % norm of the deviatoric strain
  rho_tr=2*shear.*norm_E;               % \varrho^{tr}
  p_tr=bulk.*(IOTA'*E_tr);              % trial volumetric stress
  
%
% return criteria and specification of integration points
% with plastic return to the sooth portion and to the apex of the yield
% surface
%
  denom_a= bulk.*(eta.^2);
  denom_s=shear+denom_a;
  CRIT1= rho_tr/sqrt(2) + eta.*p_tr - c ; 
  CRIT2= eta.*p_tr - denom_a.*rho_tr./(shear*sqrt(2)) - c ;
  
  % logical array indicating plastic integration points
  IND_p=CRIT1>0;
  % logical array indicating int. p. with the return to the smooth portion
  IND_s = (CRIT1>0)&(CRIT2<=0);   
  % logical array indicating integr. points with the return to the apex
  IND_a = (CRIT1>0)&(CRIT2>0);

%
% The elastic prediction of unknowns
%

  S=S_tr; DS=ELAST; 

%
% The plastic correction at the selected integration points
%

  % plastic multipliers for smooth portion of the yield surface
  lambda_s=CRIT1(IND_s)./denom_s(IND_s); 
  n_smooth=length(lambda_s);  
  % plastic multipliers for apex of the yield surface
  lambda_a=(eta(IND_a).*p_tr(IND_a)-c(IND_a))./denom_a(IND_a);
  n_apex=length(lambda_a);
      
  % correction of the stress tensors
  N_hat=dev_E(:,IND_s)./repmat(norm_E(IND_s),6,1);
  M_hat=repmat(sqrt(2)*shear(IND_s),6,1).*N_hat+IOTA*(bulk(IND_s).*eta(IND_s));
  S(:,IND_s)=S(:,IND_s)-repmat(lambda_s,6,1).*M_hat;
  S(:,IND_a)=IOTA*(c(IND_a)./eta(IND_a));

  % correction of the consistent tangent operator
  ID=DEV(:)*ones(1,n_smooth);
  NN_hat=repmat(N_hat,6,1).*kron(N_hat,ones(6,1));
  MM_hat=repmat(M_hat,6,1).*kron(M_hat,ones(6,1));
  DS(:,IND_s)=DS(:,IND_s)-...
      repmat(2*sqrt(2)*(shear(IND_s).^2).*lambda_s./rho_tr(IND_s),36,1).*(ID-NN_hat)-...
      MM_hat./repmat(denom_s(IND_s),36,1);
  DS(:,IND_a)=zeros(36,n_apex);

%     
% Update of the plastic multiplier and plastic strain (optional)
%

  if nargout<=3          
     fprintf('  plastic integration points: smooth portion=%d, apex=%d (of %d), '...
                     ,n_smooth,n_apex,n_int); 
  else 
     % plastic multiplier 
     lambda=zeros(1,n_int);  
     lambda(IND_s)=lambda_s;     % correction on the smooth portion
     lambda(IND_a)=lambda_a;     % correction at the apex  
     % plastic strain
     Ep=Ep_prev;
     Ep(:,IND_s)=Ep(:,IND_s)+([1;1;1;2;2;2]*lambda_s).*(N_hat/sqrt(2)+IOTA*(eta(IND_s)/3));
     Ep(:,IND_a)=E(:,IND_a)-IOTA*(c(IND_a)./(3*bulk(IND_a).*eta(IND_a)));     
  end

 end
