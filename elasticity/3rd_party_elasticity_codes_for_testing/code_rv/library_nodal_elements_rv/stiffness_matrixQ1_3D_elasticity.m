function [K,volumes]=stifness_matrixQ1_3D_elasticity(elements,coordinates,lambda,mu)
% Compute stiffness matrix for hexahedral mesh
%
% Code extended by tdngo

assert(size(elements,2) == 8);

%for laplace
NE=size(elements,1);        %number of elements
DIM=size(coordinates,2);    %problem dimension
 
%laplace 
NLB=8; %number of local basic functions, it must be known!
coord=zeros(DIM,NLB,NE);
for d=1:DIM
     for i=1:NLB
         coord(d,i,:)=coordinates(elements(:,i),d);
     end
end
clear coordinates

p = -1/sqrt(3);
q =  1/sqrt(3);
IP= [p, p, p;                       % Integration points. Gauss quadrature, weight = 1
     p, p, q;
     p, q, p;
     p, q, q;
     q, p, p;
     q, p, q;
     q, q, p;
     q, q, q]';

[dphi,jac] = phider(coord,IP,'Q1'); %integration rule, it must be known!  
clear coord 

jac = abs(jac);          

% Elasticity matrix
C=mu*diag([2 2 2 1 1 1]) + lambda*kron([1 0; 0 0],ones(3));

Elements=3*elements(:,kron(1:8,[1 1 1]))...
         -kron(ones(NE,1),kron(ones(1,8),[2 1 0]));
clear elements

Y=reshape(repmat(Elements,1,24)',24,24,NE);
X=permute(Y,[2 1 3]);


Z = 0;
for i = 1:size(IP,2)                      % There are 8 integration Gaussian points
    R=zeros(6,24,NE);
    R([1,4,5],1:3:22,:)=dphi(:,:,i,:);    % These indices because B = [dx 0 0; 0 dy 0; 0 0 dz; dy dx 0; dz 0 dx; 0 dz dy]
    R([4,2,6],2:3:23,:)=dphi(:,:,i,:);
    R([5,6,3],3:3:24,:)=dphi(:,:,i,:);

    Z = Z + astam(squeeze(jac(1,i,:))',amtam(R,smamt(C,permute(R,[2 1 3]))));
end
clear dphi

K=sparse(X(:),Y(:),Z(:));
volumes = sum(jac);


