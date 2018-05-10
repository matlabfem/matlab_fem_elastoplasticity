function [Xi, WF] = quadrature_volume(elem_type)

%--------------------------------------------------------------------------
% This function specifies a numerical quadrature for volume integration,
% depending on a chosen finite element. The quadratures suggested
% below can be simply replaced by another ones.
%
%  input data: 
%        elem_type - the type of Lagrange finite elements
%
%  output data:
%        Xi - local coordinates of quadrature points, size(Xi)=(2,n_q)
%        WF - weight factors, size(WF)=(1,n_q)
%--------------------------------------------------------------------------

switch(elem_type)
    
  case 'P1'      
    % - the reference triangle with coordinates:
    %               [0,0], [1,0], [0,1]
    % - 1-point quadrature rule, i.e., n_q=1
    Xi=[1/3; 1/3]; WF=1/2; 
   
  case 'P2'
    % - the reference triangle with coordinates:
    %               [0,0], [1,0], [0,1]
    % - 7-point quadrature rule, i.e., n_q=7
    Xi=[0.1012865073235, 0.7974269853531, 0.1012865073235, ...
        0.4701420641051, 0.4701420641051, 0.0597158717898, 1/3;
        0.1012865073235, 0.1012865073235, 0.7974269853531, ...
        0.0597158717898, 0.4701420641051, 0.4701420641051, 1/3];
    WF=[0.1259391805448, 0.1259391805448, 0.1259391805448, ...
        0.1323941527885, 0.1323941527885, 0.1323941527885, 0.225]/2;
        
  case 'Q1'
    % - the reference cube with coordinates:
    %         [-1,-1], [1,-1], [1,1], [-1,1]
    % - (2x2)-point quadrature rule, i.e., n_q=4 
    pt = 1/sqrt(3);
    Xi=[-pt,-pt, pt, pt 
        -pt, pt,-pt, pt];
    WF=[1,1,1,1]; 
    
  case 'Q2'
    % - the reference cube with coordinates:
    %         [-1,-1], [1,-1], [1,1], [-1,1]
    % - (3x3)-point quadrature rule, i.e., n_q=9    
    pt = sqrt(3/5);
    Xi=[-pt,  pt, pt, -pt,   0, pt,  0, -pt, 0 
        -pt, -pt, pt,  pt, -pt,  0, pt,   0, 0];
    WF=[25/81 25/81 25/81 25/81 40/81 40/81 40/81 40/81 64/81];
      
  otherwise; disp('Bad choice of element type');
end
