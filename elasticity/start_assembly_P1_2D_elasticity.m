clear level_size_P1 level_time_P1_CSV level_time_P1_RV

add_paths;

levels=8; %maximum uniform refinement level

%homogeneous material parameters
young = 206900 ;                                     % Young's modulus E
poisson =  0.29 ;                                    % Poisson's ratio nu

lambda= young*poisson/((1+poisson)*(1-2*poisson)) ;  %Lamme first parameter
mu = young./(2*(1+poisson)) ;                        %Lamme second parameter

bulk = young./(3*(1-2*poisson)) ;                    % bulk modulus K
shear = mu ;                                         % shear modulus G (equal to Lamme second parameter)    

demo=0; create_2D_mesh; %mesh for testing

for level=0:levels    
    %uniform refinement
    if (level>0)
        [coordinates,elements,dirichlet]=refinement_uniform(coordinates,elements,dirichlet);
    end
    
    figure(1); visualize_mesh;
     
    rows=numel(coordinates);
    level_size_P1(level+1)=rows; 
    
    %stiffness matrix assembly - method 1
    fprintf('technique of Cermak, Sysala and Valdman:\n')
    tic; 
    [Xi, WF] = quadrature_volume_2D('P1');    % quadrature points and weights for volume integrals
    [HatP,DHatP1,DHatP2] = local_basis_volume_2D('P1', Xi); % local basis functions and their derivatives for volume  integrals
    K=stiffness_matrix_2D(elements',coordinates',shear,bulk,DHatP1,DHatP2,WF); 
    assembly_time=toc; 
    level_time_P1_CSV(level+1)=assembly_time; 
    
    rows=size(K,1);
    fprintf('level=%d, ', level);
    fprintf('time spent on K: %6.1e seconds, ',assembly_time(end));
    fprintf('rows of matrix =%d ',rows);
    fprintf('\n');  
    
    %stiffness matrix assembly - method 2
    fprintf('technique of Rahman and Valdman: \n')
    tic; 
    K=stiffness_matrixP1_2D_elasticity(elements,coordinates,lambda,mu); 
    assembly_time=toc; 
    level_time_P1_RV(level+1)=assembly_time; 
    
    fprintf('level=%d, ', level);
    fprintf('time spent on K: %6.1e seconds, ',assembly_time(end));
    fprintf('rows of matrix =%d ',rows);
    fprintf('\n');  
    
    fprintf('-----------------------------------------------\n')
end

% output information
fprintf('\n')
for level=0:levels
    fprintf('%d ', level);
    fprintf('& ');
    fprintf('%d ', level_size_P1(level+1));
    fprintf('& ');
    fprintf('%2.2f ', level_time_P1_CSV(level+1));
    fprintf('& ');
    fprintf('%2.2f ', level_time_P1_RV(level+1));
    fprintf('\\\\');
    fprintf('\n');
end


