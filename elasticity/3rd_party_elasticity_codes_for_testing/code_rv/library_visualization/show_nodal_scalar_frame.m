function h=show_nodal_scalar_frame(nodalValue,elements,coordinates,nodalDisplacement)
switch nargin, 
    case 3,
        nodalDisplacement=0*coordinates;
    case {0, 1, 2}
        fprintf('missing parameters')
end

if size(coordinates,2)==2
    X=coordinates(:,1)+nodalDisplacement(:,1);
    Y=coordinates(:,2)+nodalDisplacement(:,2);
    
    h=trisurf(elements,X,Y,nodalValue,'FaceColor','interp');
    
    %XX=X(elements)';
    %YY=Y(elements)';
    %ZZ=nodalValue(elements)';
    %CC=ZZ;
    %fill3(XX,YY,ZZ,CC,'FaceColor','interp','Linewidth',0.1);
elseif size(coordinates,2)==3  %3D only testing
    X=coordinates(:,1)+nodalDisplacement(:,1);
    Y=coordinates(:,2)+nodalDisplacement(:,2);
    Z=coordinates(:,3)+nodalDisplacement(:,3);
    
    elements=[elements(:,[1 2 3]); elements(:,[1 2 4]); elements(:,[2 3 4]); elements(:,[3 1 4])];
    
    h=trisurf(elements,X,Y,Z,nodalValue,'FaceColor','interp');
    
    %XX=X(elements)';
    %YY=Y(elements)';
    %ZZ=Z(elements)';
    %CC=nodalValue(elements)';
    %fill3(XX,YY,ZZ,CC,'FaceColor','interp','Linewidth',0.1);
end

%set(gcf,'renderer','zbuffer');