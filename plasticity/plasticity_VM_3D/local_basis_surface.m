function [HatP_s,DHatP1_s,DHatP2_s] = local_basis_surface(elem_type, Xi_s)

%--------------------------------------------------------------------------
% This function evaluates local basis functions on a surface and their 
% derivatives at prescribed quadrature points depending on a chosen
% finite elements.
%
%  input data: 
%       elem_type - the type of Lagrange finite elements
%       Xi_s      - coordinates of the quadrature points, size(Xi_s)=(3,n_q_s)
%
%  output data:
%       HatP_s   - values of basis functions at the quadrature points,
%                  size(HatP_s)=(n_p_s,n_q_s)
%       DHatP1_s - derivatives of basis functions at the quadrature points 
%                  in the direction xi_1, size(DHatP1_s)=(n_p_s,n_q_s)
%       DHatP2_s - derivatives of basis functions at the quadrature points 
%                  in the direction xi_2, size(DHatP2_s)=(n_p_s,n_q_s)
%       n_p_s    - number of basis functions on the surface
%       n_q_s    - number of integration points within the surface
%--------------------------------------------------------------------------

xi_1 = Xi_s(1,:); xi_2 = Xi_s(2,:); 

switch(elem_type)
  case 'P1'
    % - the reference triangle with coordinates:
    %               [0,0], [1,0], [0,1]
    % - n_p_s=3, n_q_s=length(xi_1)
    HatP_s = [1-xi_1-xi_2; xi_1; xi_2] ;
    DHatP1_s = [ -1; 1; 0];
    DHatP2_s = [ -1; 0; 1];
    
  case 'P2'
    % - the reference triangle with coordinates:
    %               [0,0], [1,0], [0,1]
    % - coordinates of midpoints:
    %             [1/2,1/2], [0,1/2], [1/2,0], 
    % - n_p_s=6, n_q_s=length(xi_1)
    % - barycentric coordinates: xi_0, xi_1, xi_2, where
    xi_0 = 1-xi_1-xi_2;
    n_q_s = length(xi_1);
    HatP_s = [ xi_0.*(2*xi_0-1); xi_1.*(2*xi_1-1); xi_2.*(2*xi_2-1);
                   4*xi_1.*xi_2;     4*xi_0.*xi_2;     4*xi_0.*xi_1 ];
    DHatP1_s = [ -4*xi_0+1; 4*xi_1-1; zeros(1,n_q_s); 
                    4*xi_2;  -4*xi_2;  4*(xi_0-xi_1) ];
    DHatP2_s = [ -4*xi_0+1; zeros(1,n_q_s); 4*xi_2-1; 
                    4*xi_1;  4*(xi_0-xi_2);  -4*xi_1 ];
       
  case 'Q1'
    % - the reference square with coordinates:
    %         [-1,-1], [1,-1], [1,1], [-1,1] 
    % - n_p_s=4, n_q_s=length(xi_1)    
    HatP_s=[(1-xi_1).*(1-xi_2)/4; (1+xi_1).*(1-xi_2)/4;
            (1+xi_1).*(1+xi_2)/4; (1-xi_1).*(1+xi_2)/4 ];
    DHatP1_s = [-(1-xi_2)/4;  (1-xi_2)/4;
                 (1+xi_2)/4; -(1+xi_2)/4 ];
    DHatP2_s = [-(1-xi_1)/4; -(1+xi_1)/4;
                 (1+xi_1)/4;  (1-xi_1)/4 ];
     
  case 'Q2'
    % - the reference square with coordinates:
    %         [-1,-1], [1,-1], [1,1], [-1,1]
    % - coordinates of midpoints:
    %          [0,-1],  [1,0], [0,1], [-1,0]
    % - n_p_s=8, n_q_s=length(xi_1)    
    HatP_s = [ (1-xi_1).*(1-xi_2).*(-1-xi_1-xi_2)/4;
               (1+xi_1).*(1-xi_2).*(-1+xi_1-xi_2)/4;
               (1+xi_1).*(1+xi_2).*(-1+xi_1+xi_2)/4;
               (1-xi_1).*(1+xi_2).*(-1-xi_1+xi_2)/4;            
                            (1-xi_1.^2).*(1-xi_2)/2;
                            (1+xi_1).*(1-xi_2.^2)/2;
                            (1-xi_1.^2).*(1+xi_2)/2;
                            (1-xi_1).*(1-xi_2.^2)/2 ];

    DHatP1_s = [ (1-xi_2).*(2*xi_1+xi_2)/4;
                 (1-xi_2).*(2*xi_1-xi_2)/4;
                 (1+xi_2).*(2*xi_1+xi_2)/4;
                 (1+xi_2).*(2*xi_1-xi_2)/4;
                           -xi_1.*(1-xi_2);
                             (1-xi_2.^2)/2;
                           -xi_1.*(1+xi_2);
                            -(1-xi_2.^2)/2 ];

    DHatP2_s = [ (1-xi_1).*( xi_1+2*xi_2)/4;
                 (1+xi_1).*(-xi_1+2*xi_2)/4;
                 (1+xi_1).*( xi_1+2*xi_2)/4;
                 (1-xi_1).*(-xi_1+2*xi_2)/4;
                             -(1-xi_1.^2)/2;
                            -(1+xi_1).*xi_2;
                              (1-xi_1.^2)/2;
                            -(1-xi_1).*xi_2 ];
    
  otherwise; disp('Bad choise of element type');
end
