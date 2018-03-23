function show_mesh(elements,coordinates)
if size(coordinates,2)==2    %2D objects 
   %works for both triangles and rectangles         
        X=reshape(coordinates(elements',1),size(elements,2),size(elements,1));
        Y=reshape(coordinates(elements',2),size(elements,2),size(elements,1));        
        fill(X,Y,'w'); 
elseif size(coordinates,2)==3   %3D object
    if size(elements,2)==4      %tetrahedra
        tetramesh(elements,coordinates,'FaceAlpha',1,'FaceColor','white');%camorbit(20,0);
    end
    
    if size(elements,2)==8      %hexahedra
       faces1=elements(:,1:4);
       faces2=elements(:,5:8);
       faces3=elements(:,[1 4 8 5]);
       faces4=elements(:,[2 3 7 6]);
       faces5=elements(:,[1 2 6 5]);
       faces6=elements(:,[3 4 8 7]);
       
       faces=[faces1; faces2; faces3; faces4; faces5; faces6];
       
       X=reshape(coordinates(faces',1),size(faces,2),size(faces,1));
       Y=reshape(coordinates(faces',2),size(faces,2),size(faces,1)); 
       Z=reshape(coordinates(faces',3),size(faces,2),size(faces,1)); 
       
       fill3(X,Y,Z,[0.0 0.0 0.0],'FaceAlpha',1,'FaceColor','white');     
    end
    
    if size(elements,2)==3      % triangular faces
       faces=elements;
       
       X=reshape(coordinates(faces',1),size(faces,2),size(faces,1));
       Y=reshape(coordinates(faces',2),size(faces,2),size(faces,1)); 
       Z=reshape(coordinates(faces',3),size(faces,2),size(faces,1)); 
       
       fill3(X,Y,Z,[0.3 0.3 0.9]);     
    end
    
end


