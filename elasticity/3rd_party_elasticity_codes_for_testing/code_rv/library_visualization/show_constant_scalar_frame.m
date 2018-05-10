function show_constant_scalar_frame(constantValue,elements,coordinates,nodalDisplacement)
    
if (nargin==4)     
    if norm(nodalDisplacement)==0
       factor=1;
    else
       factor =10^(-round(log10(max(max(nodalDisplacement)))));
    end   
else
    nodalDisplacement=zeros(size(coordinates));
end

if size(coordinates,2)==2
    X=coordinates(:,1)+nodalDisplacement(:,1);
    Y=coordinates(:,2)+nodalDisplacement(:,2);
    
    XX=X(elements)';
    YY=Y(elements)';
    ZZ=kron(ones(1,size(elements,2)),constantValue)';
    CC=ZZ;
    fill3(XX,YY,ZZ,CC,'FaceColor','interp','Linewidth',0.1);

if size(coordinates,2)==3  %3D only testing
    
    X=coordinates(:,1)+nodalDisplacement(:,1);
    Y=coordinates(:,2)+nodalDisplacement(:,2);
    Z=coordinates(:,3)+nodalDisplacement(:,3);
    
    elements=[elements(:,[1 2 3]); elements(:,[1 2 4]); elements(:,[2 3 4]); elements(:,[3 1 4])];
    constantValue=[constantValue; constantValue; constantValue; constantValue];
    
    XX=X(elements)';
    YY=Y(elements)';
    ZZ=Z(elements)';
    CC=constantValue(elements)';
    fill3(XX,YY,ZZ,CC,'FaceColor','interp','Linewidth',0.1);
    
    trimesh(elements,X,Y,Z);
end



end
