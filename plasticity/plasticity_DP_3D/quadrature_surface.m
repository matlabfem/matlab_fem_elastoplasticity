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
%       Xi_s - local coordinates of quadrature points, size(Xi_s)=(2,n_q_s)
%       WF_s - weight factors, size(WF_s)=(1,n_q_s)
%--------------------------------------------------------------------------

switch(elem_type)
    
  case 'P1'      
    % - surface is created by the reference triangle with coordinates:
    %               [0,0], [1,0], [0,1]
    % - 1-point quadrature rule, i.e., n_q_s=1
    Xi_s=[1/3; 1/3]; WF_s=1/2; 
   
  case 'P2'
    % - surface is created by the reference triangle with coordinates:
    %               [0,0], [1,0], [0,1]
    % - 3-point quadrature rule, i.e., n_q_s=3
    Xi_s=[1/6, 1/6, 2/3
          1/6, 2/3, 1/6];
    WF_s=[1/6, 1/6, 1/6]; 
        
  case 'Q1'
    % - surface is created by the reference square with coordinates:
    %         [-1,-1], [1,-1], [1,1], [-1,1]
    % - (2x2)-point quadrature rule, i.e., n_q_s=4  
    pt = 1/sqrt(3);
    Xi_s=[-pt,-pt, pt, pt 
          -pt, pt,-pt, pt];
    WF_s=[1,1,1,1]; 
    
  case 'Q2'
    % - surface is created by the reference square with coordinates:
    %         [-1,-1], [1,-1], [1,1], [-1,1]
    % - (2x2)-point quadrature rule, i.e., n_q_s=4  
    pt = 1/sqrt(3);
    Xi_s=[-pt,-pt, pt, pt 
          -pt, pt,-pt, pt];
    WF_s=[1,1,1,1]; 
      
  otherwise; disp('Bad choise of element type');
end
