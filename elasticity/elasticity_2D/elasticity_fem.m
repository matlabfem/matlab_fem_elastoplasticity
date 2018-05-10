
function [assembly_time,rows]=elasticity_fem(elem_type,level,draw)

% =========================================================================
%
%  This function constructs the mesh, assembles the elastic stiffness
%  matrix and the load vector, and computes displacements depending on 
%  an element type and mesh density. Optionally, it enables to visualize 
%  the results and the mesh.
%
%  input data:
%    elem_type - the type of finite elements; available choices:
%                'P1', 'P2', 'Q1', 'Q2'
%    level     - an integer defining mesh density
%    draw      - a logical value enbling visualization of the results
%
%  output data:
%    assembly_time - assembly time for the elastic stiffness matrix
%    rows          - number of rows in the stiffness matrix
%
% ======================================================================
%

%
% specification of input data to the function
%
  if nargin==0
    elem_type='P1';
  end
  if nargin<=1
     level=1;
  end
  if nargin<=2
    draw = 1;
  end
  
%
% The main input data 
%
  
  % values of elastic material parameters                     
  young = 206900;          % Young's modulus
  poisson =  0.29;         % Poisson's ratio
  shear = young/(2*(1+poisson)) ;        % shear modulus
  bulk = young/(3*(1-2*poisson)) ;       % bulk modulus
  
  % constant volume forces in each directions 
  volume_force = [0,-1] ; 
  
  % constant tranction on the back side of the body in each direction
  traction_force = [0,450] ; 

%
% mesh generation
%

  % geometrical parameters (choose only integers, size_hole < size_xy)
  size_xy = 10;      % size of the body in direction x and y 
  size_hole = 5;     % size of the hole in the body
                     % the projection of the hole to the xy-plane is square 
                     
  % the mesh generation depending prescribed finite elements        
  switch(elem_type)
    case 'P1'; [COORD,ELEM,SURF,NEUMANN,DIRICHLET,Q]=mesh_P1(level,size_xy,size_hole);
    case 'P2'; [COORD,ELEM,SURF,NEUMANN,DIRICHLET,Q]=mesh_P2(level,size_xy,size_hole);
    case 'Q1'; [COORD,ELEM,SURF,NEUMANN,DIRICHLET,Q]=mesh_Q1(level,size_xy,size_hole);
    case 'Q2'; [COORD,ELEM,SURF,NEUMANN,DIRICHLET,Q]=mesh_Q2(level,size_xy,size_hole);
    otherwise
          disp('bad choice of element type');
  end  

%
% Data from the reference element
%
  
  % quadrature points and weights for volume and surface integration
  [Xi, WF] = quadrature_volume(elem_type);
  [Xi_s, WF_s] = quadrature_surface(elem_type);
  
  % local basis functions and their derivatives for volume and surface
  [HatP,DHatP1,DHatP2] = local_basis_volume(elem_type, Xi);
  [HatP_s,DHatP1_s] = local_basis_surface(elem_type, Xi_s);

% Assembling of the stiffness matrix
  
  % auxiliary values
  n_e=size(ELEM,2);     % number of elements
  n_q=length(WF);       % number of quadratic points
  n_int = n_e*n_q ;     % total number of integrations points  
 
  % values of material parameters at integration points      
  shear =shear*ones(1,n_int);
  bulk=bulk*ones(1,n_int);

  % stiffness matrix assembly and the assembly time
  tic; 
  [K,WEIGHT]=elastic_stiffness_matrix(ELEM,COORD,shear,bulk,DHatP1,DHatP2,WF);  
  assembly_time=toc; 
  rows=size(K,1);
  fprintf('level=%d, ', level);
  fprintf('time spent on K: %6.1e seconds, ',assembly_time(end));
  fprintf('rows of matrix =%d ',rows);
  fprintf('\n');  
  
%
% Assembling of the vector of volume forces
%     

  % values of the density function f_V at integration points
  f_V_int=volume_force'*ones(1,n_int); % size(f_V_int)=(2,n_int)
  
  % assembling of the vector of volume forces
  f_V=vector_volume(ELEM,COORD,f_V_int,HatP,WEIGHT);
  
%
% Assembling of the vector of traction forces
%  
  
  % values of the density function f_t at surface integration points
  n_e_s=size(NEUMANN,2); % number of surface elements
  n_q_s=length(WF_s);   % number of quadrature points on a surface element
  n_int_s=n_e_s*n_q_s ; % number of integration points on the surface
                        % (on the upper side of the body)
  f_t_int=traction_force'*ones(1,n_int_s); % size(f_V_int)=(2,n_int_s)
  
  % assembling of the vector of traction (surface) forces
  f_t=vector_traction(NEUMANN,COORD,f_t_int,HatP_s,DHatP1_s,WF_s);
  
%
% Ud - 2*n_n array representing constant displacement in y-direction 
% on the front of the body.
%
  Ud = 0.5*DIRICHLET;
  
%
% Processing
%    
  
  % the right hand side in a system of linear equations
  f = f_t(:) + f_V(:) - K*Ud(:);
  
  % computation of displacements
  U = Ud;
  U(Q) = K(Q,Q)\f(Q) ;
  
  % stored energy
  e=(1/2)*U(:)'*K*U(:)-(f_t(:)'+f_V(:)')*U(:); 
  fprintf('stored energy = %d \n',e);

%  
% visualization of the mesh and displacements
%
  if draw
      
    % mesh visualization
    draw_mesh(COORD,ELEM,elem_type)      
     
    % total displacements + deformed shape
    U_total = sqrt(U(1,:).^2 + U(2,:).^2);
    draw_displacement(COORD,ELEM,U,U_total,elem_type,size_xy,size_hole)
      
  end

end
