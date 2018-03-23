function [A,areas]=stiffness_matrixP1_2D_elasticity(elements,coordinates,lambda,mu)

%for laplace
NE=size(elements,1); %number of elements
DIM=size(coordinates,2); %problem dimension
 
%laplace 
NLB=3; %number of local basic functions, it must be known!
coord=zeros(DIM,NLB,NE);
for d=1:DIM
     for i=1:NLB
         coord(d,i,:)=coordinates(elements(:,i),d);
     end
end   
IP=[1/3;1/3];
[dphi,jac] = phider(coord,IP,'P1'); %integration rule, it must be known! 
dphi = squeeze(dphi); 
areas=abs(squeeze(jac))/factorial(DIM);

%elasticity matrix is derived from laplace matrix
R=zeros(3,6,NE);
R([1,3],[1,3,5],:)=dphi;  
R([3,2],[2,4,6],:)=dphi;
clear dphi

C=mu*[2 0 0;0 2 0;0 0 1] + lambda*[1 1 0;1 1 0;0 0 0];

Elements=2*elements(:,[1 1 2 2 3 3])-kron(ones(NE,1),[1,0,1,0,1,0]);
NLB=6;
Y=reshape(repmat(Elements,1,NLB)',NLB,NLB,NE);
X=permute(Y,[2 1 3]);
Z=astam(areas',amtam(R,smamt(C,permute(R,[2 1 3]))));

A=sparse(X(:),Y(:),Z(:));

end


