function [Xi_s, WF_s] = quadrature_surface(elem_type)

%--------------------------------------------------------------------------
% This function specifies a numerical quadrature for surface integration,
% depending on a chosen finite element. The quadratures suggested
% below can be simply replaced by another ones.
%
%  input data: 
%       elem_type - the type of Lagrange finite elements
%
%  output data:
%       Xi_s - local coordinates of quadrature points, size(Xi_s)=(1,n_q_s)
%       WF_s - weight factors, size(WF_s)=(1,n_q_s)
%--------------------------------------------------------------------------

switch(elem_type)
    
  case {'P1','Q1'}      
    % - surface is created by the reference line with coordinates:
    %               -1, 1
    % - 1-point quadrature rule, i.e., n_q_s=1
    Xi_s=0; WF_s=2; 
   
  case {'P2','Q2'}
    % - surface is created by the reference line with coordinates:
    %               -1, 1
    % - 2-point quadrature rule, i.e., n_q_s=2
    pt = 1/sqrt(3);
    Xi_s=[-pt,pt];
    WF_s=[1,1]; 
      
  otherwise; disp('Bad choise of element type');
end
