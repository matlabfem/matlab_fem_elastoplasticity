function [HatP_s,DHatP1_s] = local_basis_surface(elem_type, Xi_s)

%--------------------------------------------------------------------------
% This function evaluates local basis functions on a surface and their 
% derivatives at prescribed quadrature points depending on a chosen
% finite elements.
%
%  input data: 
%       elem_type - the type of Lagrange finite elements
%       Xi_s - coordinates of the quadrature points, size(Xi_s)=(1,n_q_s)
%
%  output data:
%       HatP_s - values of basis functions at the quadrature points,
%              size(HatP_s)=(n_p_s,n_q_s)
%       DHatP1_s - derivatives of basis functions at the quadrature points 
%                in the direction xi_1, size(DHatP1_s)=(n_p_s,n_q_s)
%       n_p_s - number of basis functions on the surface
%       n_q_s - number of integration points within the surface
%--------------------------------------------------------------------------

xi = Xi_s(1,:); 

switch(elem_type)
  case {'P1','Q1'}
    % - the reference line with coordinates:
    %               -1, 1
    % - n_p_s=2, n_q_s=length(xi_1)
    HatP_s = (1/2)*[1-xi; 1+xi] ;
    DHatP1_s = [-1/2; 1/2];
    
  case {'P2','Q2'}
    % - the reference triangle with coordinates:
    %               -1, 1
    % - coordinates of midpoint:
    %             0
    % - n_p_s=3, n_q_s=length(xi_1)
    HatP_s = [xi.*(xi-1)/2; xi.*(xi+1)/2; (xi+1).*(1-xi)] ;
    DHatP1_s = [ xi-1/2; xi+1/2; -2*xi];
       
  otherwise; disp('Bad choise of element type');
end
