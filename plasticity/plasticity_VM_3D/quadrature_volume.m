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
%        Xi - local coordinates of quadrature points, size(Xi)=(3,n_q)
%        WF - weight factors, size(WF)=(1,n_q)
%--------------------------------------------------------------------------

switch(elem_type)
    
  case 'P1'      
    % - the reference tetrahedron with coordinates:
    %               [0,0,0], [1,0,0], [0,1,0], [0,0,1]
    % - 1-point quadrature rule, i.e., n_q=1
    Xi=[1/4; 1/4; 1/4]; WF=1/6; 
   
  case 'P2'
    % - the reference tetrahedron with coordinates:
    %               [0,0,0], [1,0,0], [0,1,0], [0,0,1]
    % - 11-point quadrature rule, i.e., n_q=11
    Xi=[1/4, 0.0714285714285714, 0.785714285714286,  0.0714285714285714, 0.0714285714285714, 0.399403576166799, 0.100596423833201, 0.100596423833201, 0.399403576166799, 0.399403576166799, 0.100596423833201
        1/4, 0.0714285714285714, 0.0714285714285714, 0.785714285714286,  0.0714285714285714, 0.100596423833201, 0.399403576166799, 0.100596423833201, 0.399403576166799, 0.100596423833201, 0.399403576166799
        1/4, 0.0714285714285714, 0.0714285714285714, 0.0714285714285714, 0.785714285714286,  0.100596423833201, 0.100596423833201, 0.399403576166799, 0.100596423833201, 0.399403576166799, 0.399403576166799];
    WF=[-0.013155555555555, 0.007622222222222*ones(1,4), 0.024888888888888*ones(1,6)]; 
        
  case 'Q1'
    % - the reference cube with coordinates:
    %         [-1,-1,-1], [1,-1,-1], [1,1,-1], [-1,1,-1], 
    %          [-1,-1,1],  [1,-1,1],  [1,1,1],  [-1,1,1]
    % - (2x2x2)-point quadrature rule, i.e., n_q=8  
    pt = 1/sqrt(3);
    Xi=[-pt,-pt,-pt,-pt, pt, pt, pt, pt
        -pt,-pt, pt, pt,-pt,-pt, pt, pt
        -pt, pt,-pt, pt,-pt, pt,-pt, pt];
    WF=[1,1,1,1,1,1,1,1]; 
    
  case 'Q2'
    % - the reference cube with coordinates:
    %         [-1,-1,-1], [1,-1,-1], [1,1,-1], [-1,1,-1], 
    %          [-1,-1,1],  [1,-1,1],  [1,1,1],  [-1,1,1]
    % - (3x3x3)-point quadrature rule, i.e., n_q=27    
    pt = sqrt(3/5);
    Xi=[-pt,-pt,-pt,-pt,-pt,-pt,-pt,-pt,-pt,  0,  0,  0,  0,  0,  0,  0,  0,  0, pt, pt, pt, pt, pt, pt, pt, pt, pt
        -pt,-pt,-pt,  0,  0,  0, pt, pt, pt,-pt,-pt,-pt,  0,  0,  0, pt, pt, pt,-pt,-pt,-pt,  0,  0,  0, pt, pt, pt
        -pt,  0, pt,-pt,  0, pt,-pt,  0, pt,-pt,  0, pt,-pt,  0, pt,-pt,  0, pt,-pt,  0, pt,-pt,  0, pt,-pt,  0, pt];
    WF=    ((5/9)^3)* [1 0 1 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 1 0 1] +...
     ((5/9)^2)*(8/9)* [0 1 0 1 0 1 0 1 0 1 0 1 0 0 0 1 0 1 0 1 0 1 0 1 0 1 0] +...     
     ((8/9)^2)*(5/9)* [0 0 0 0 1 0 0 0 0 0 1 0 1 0 1 0 1 0 0 0 0 0 1 0 0 0 0] +...
           ((8/9)^3)* [0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0];
      
  otherwise; disp('Bad choice of element type');
end
