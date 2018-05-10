function [A,areas,Z,X,Y]=stiffness_matrixP2_2D(elements,coordinates,coeffs)
%coeffs can be either P0 (elementwise constant) or P1 (elementwise nodal) function 
%represented by a collumn vector with size(elements,1) or size(coordinates,1) entries
%if coeffs is not provided then coeffs=1 is assumed globally

NE=size(elements,1); %number of elements
DIM=size(coordinates,2); %problem dimension

%particular part for a given element in a given dimension
NLB=6; %number of local basic functions, it must be known!
coord=zeros(DIM,NLB,NE);
for d=1:DIM
    for i=1:NLB
        coord(d,i,:)=coordinates(elements(:,i),d);
    end
end  

p=3; [IP,w]=inttri(p); IP=IP'; w=w*2;

[dphi,jac] = phider(coord,IP,'P2'); 

areas=abs(squeeze(jac))/factorial(DIM);

if p>=2
   areas=areas(1,:)';
end
    

if (nargin<3)
    product=zeros(NLB,NLB,NE);
    for i=1:size(IP,2)
        product=product+w(i)*amtam(squeeze(dphi(:,:,i,:)),squeeze(dphi(:,:,i,:)));
    end
    Z=astam(areas,product);  
else
    if numel(coeffs)==size(coordinates,1)  %P1->P0 averaging
        coeffs=evaluate_average_point(elements,coeffs);
    end  
    Z=astam((areas.*coeffs)',amtam( dphi,dphi));
end
Y=reshape(repmat(elements,1,NLB)',NLB,NLB,NE);

%copy this part for a creation of a new element
X=permute(Y,[2 1 3]);
A=sparse(X(:),Y(:),Z(:));  