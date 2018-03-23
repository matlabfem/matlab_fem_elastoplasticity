%
% Elasticité linéaire 2D: Membrane de Cook
% Solution par élimination de Gauss
%----------------------------------------------------------
% constantes
E=30000; nu=0.4;
Penalty=10^15*E;

% chargement du maillage
load cook995 p t ibcd ibcneum
n=size(p,1); nn=2*n;

% d.d.l. de Dirichlet 
nnb=2*length(ibcd);
ibc=zeros(nnb,1); ibc(1:2:end)=2*ibcd-1; ibc(2:2:end)=2*ibcd;

% assemblage de la matrice et second membre
K=kelas2drgd(p,t,E,nu);
f1=zeros(n,1); f2=-.75*ones(n,1);
f=kelas2drhs(p,t,f1,f2);
g=kelas2drhsn(p,ibcneum,0,10);
clear ibcn ibcn1 ibcn2

% conditions aux limites de Dirichlet par pénalisation
K(ibc,ibc)=K(ibc,ibc)+Penalty*speye(nnb);
b=f+g; b(ibc)=0;

% solution par élimination de Gauss
u=K\b;

% contraintes de Von Mises, visualisation 
s=kelas2dvmes(p,t,u,E,nu);
colormap(1-gray), kelas2dshow(p,t,u,20,s)
