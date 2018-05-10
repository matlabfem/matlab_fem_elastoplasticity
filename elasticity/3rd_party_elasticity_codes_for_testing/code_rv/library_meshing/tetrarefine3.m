function [XYZ,Elm,ERef] = tetrarefine3(XYZ,Elm,ERef);
% author: Ales Janka, ales.janka@unifr.ch, http://perso.unifr.ch/ales.janka
% function [XYZ,Elm,ERef] = tetrarefine3(XYZ,Elm,ERef);
% refines tetrahedras bu cutting each edge in half and making 8 new 
% finer tetrahedra out of one corser one. 
% old nodal coords come first in XYZ, then the new ones
% The new tetrahedra are similar to the old one, 
% no degeneration is supposed to occur as at most 3 congruence classes of 
% tetrahedra appear, even when re-applied iteratively (provided that 
% Elm() is not modified between two applications - ordering of vertices 
% in tetrahedra matters not only for positivity of volumes).
%*** References:
% Juergen Bey: Simplicial grid refinement: on Freudenthal s algorithm and 
%   the optimal number of congruence classes, Numer.Math. 85 (2000), 
%   no. 1, 1--29, or
% Juergen Bey: Tetrahedral grid refinement, Computing 55 (1995), 
%   no. 4, 355--378, or
% http://citeseer.ist.psu.edu/bey95tetrahedral.html
%*** Obsoletes the observations in:
% Carre G., Carte G., Guillard H., Lanteri S.: Multigrid strategies for 
%   CFD problems on non-structured meshes, Multigrid methods, VI 
%   (European MG conference Gent 1999), 1--10, 
%   Lect.Notes.Comput.Sci.Engrg. 190 (2000), no. 11-12, 1467--1482.

nelm = size(Elm,1);
nnod = size(XYZ,1);

if (nargin<3)
  ERef=[];
end;

%fprintf('--- midsides begin\n');
edge = [1 2;
        1 3;
        1 4;
        2 3;
        2 4;
        3 4];
Edges = Elm(:,edge(1,:));
[Edges,I,J] = unique(sort(Edges,2),'rows');
clear I
for i=2:6
  n = size(Edges,1);
  Edges = [Edges;Elm(:,edge(i,:))];
  n1 = size(Edges,1);
  [Edges,I1,J1] = unique(sort(Edges,2),'rows');
  clear I1
  J = [J1(J);J1(n+1:n1,1)];
  clear J1
end;
%fprintf('--- midsides done\n');

XYZ  = [XYZ; (XYZ(Edges(:,1),:) + XYZ(Edges(:,2),:))/2];
J = J+nnod;

clear Edges

Elm = [Elm, reshape(J,nelm,6)];
clear J

%fprintf('--- elements begin\n');
Elm = [Elm(:, 1) Elm(:, 5) Elm(:, 6) Elm(:, 7), ...
       Elm(:, 5) Elm(:, 2) Elm(:, 8) Elm(:, 9), ...
       Elm(:, 6) Elm(:, 8) Elm(:, 3) Elm(:, 10), ...
       Elm(:, 7) Elm(:, 9) Elm(:,10) Elm(:, 4), ...
       Elm(:, 5) Elm(:, 6) Elm(:, 7) Elm(:, 9), ...
       Elm(:, 8) Elm(:, 6) Elm(:, 5) Elm(:, 9), ...
       Elm(:, 6) Elm(:, 7) Elm(:, 9) Elm(:,10), ...
       Elm(:, 9) Elm(:, 8) Elm(:, 6) Elm(:,10)];
Elm = Elm';
Elm = reshape(Elm,4,8*nelm);
Elm = Elm';
if ~isempty(ERef);
  ERef = ERef * [1 1 1 1 1 1 1 1];
  ERef = reshape(ERef',8*nelm,1);
end;
%fprintf('--- elements done\n');

