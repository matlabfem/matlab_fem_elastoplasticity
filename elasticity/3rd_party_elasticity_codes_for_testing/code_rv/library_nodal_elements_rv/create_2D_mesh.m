if demo==0
   %for demo=0
    %coarse mesh of the unit square (-0.5, 0.5)^2 (with half middle edge boundary)
    coordinates=0.5*[0 0; 1 0; 2 0; 2 1; 1 1; 0 1; 0 2; 1 2; 2 2];
    elements=[1 2 5; 5 6 1; 2 3 4; 2 4 5; 6 5 8; 6 8 7; 5 4 9; 5 9 8]; 
    dirichlet=[1 2; 2 3; 3 4; 4 5; 4 9; 9 8; 8 7; 7 6; 6 1];

    %coarse mesh of the L-shape
    coordinates=[0 0; 1 0; 2 0; 2 1; 1 1; 0 1; 0 2; 1 2];
    elements=[1 2 5; 5 6 1; 2 3 4; 2 4 5; 6 5 8; 6 8 7];   

    %coarse mesh of the L-shape (rotated)
    coordinates=[0 -2; -1 -1; 1 -1; 0 0; 2 0; 1 1; -1 1; 0 2];
    elements=[1 3 2; 3 4 2; 3 5 4; 5 6 4; 4 6 7; 6 8 7];
    dirichlet=[1 2; 2 3; 3 4; 4 5; 5 8; 8 7; 7 6; 6 1];
    neumann=[]; 
else   
    %coarse mesh of the unit square (0,1)^2
    coordinates=[0 0; 1 0; 0 1; 1 1];
    elements=[1 2 3; 2 4 3];   
    %dirichlet=[1 2; 2 4; 4 3; 3 1]; 

    if (demo==1)
        dirichlet=[1 2; 2 4; 4 3; 3 1];
        neumann=[];
    else  %for demo>1
        dirichlet=[2 4; 4 3; 3 1]; 
        neumann=[1 2];
    end  
end

