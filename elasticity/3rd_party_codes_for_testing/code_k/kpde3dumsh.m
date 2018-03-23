function [p,t,bpx,bpy,bpz]=kpde3dumsh(ax,bx,ay,by,az,bz,nx,ny,nz)
%KPDE3DUMSH Uniform mesh generation of (ax,bx)*(ay,by)*(az,bz)
%--------------------------------------------------------------------
%      [p,t,bpx,bpy,bpz]=kpde3dumsh(ax,bx,ay,by,az,bz,nx,ny,nz)
% 
% Input:
%     (ax,bx) : x-interval
%     (ay,by) : y-interval
%     (az,bz) : z-interval
%    nx,ny,nz : number of points in x, y and z
% Ouput:
%          p  : Nodes coordinates, array np*3
%          t  : Tetrahedron vertices, array nt*4
% bpx,bpy,bpz : boundary nodes, arrays nx*2, ny*2 and nz*2
%               bpx(:,1), bpx(:,2) boundary nodes at x=ax and x=bx  
%               bpy(:,1), bpy(:,2) boundary nodes at y=ay and y=by  
%               bpz(:,1), bpz(:,2) boundary nodes at z=az and z=bz  
%--------------------------------------------------------------------
np=nx*ny*nz; nq=(nx-1)*(ny-1)*(nz-1); nt=6*nq;                   

% nodes coordinates
hx=(bx-ax)/(nx-1); hy=(by-ay)/(ny-1); hz=(bz-az)/(nz-1);
[x,y,z]=ndgrid(ax:hx:bx,ay:hy:by,az:hz:bz);
p=[x(:),y(:),z(:)];

% quadrangles
ip=[1:nx*ny*nz]'; ijk=reshape(ip,[nx,ny,nz]);
b1=squeeze(ijk(:,1,:));  ib1=b1(:);
b2=squeeze(ijk(nx,:,:)); ib2=b2(:);
b3=squeeze(ijk(:,ny,:)); ib3=b3(:);
b4=squeeze(ijk(1,:,:));  ib4=b4(:);
b5=squeeze(ijk(:,:,1));  ib5=b5(:);
b6=squeeze(ijk(:,:,nz)); ib6=b6(:);

ib236=union(ib2,union(ib3,ib6)); ip1=setdiff(ip,ib236);
ib346=union(ib3,union(ib4,ib6)); ip2=setdiff(ip,ib346);
ib345=union(ib3,union(ib4,ib5)); ip3=setdiff(ip,ib345);
ib235=union(ib2,union(ib3,ib5)); ip4=setdiff(ip,ib235);
ib126=union(ib1,union(ib2,ib6)); ip5=setdiff(ip,ib126);
ib146=union(ib1,union(ib4,ib6)); ip6=setdiff(ip,ib146);
ib145=union(ib1,union(ib4,ib5)); ip7=setdiff(ip,ib145);
ib125=union(ib1,union(ib2,ib5)); ip8=setdiff(ip,ib125);

% triangulation
t=zeros(nt,4);
iq1=1:nq; iq2=iq1+nq; iq3=iq2+nq;  
iq4=iq3+nq; iq5=iq4+nq; iq6=iq5+nq; 
t(iq1,1)=ip1; t(iq1,2)=ip2; t(iq1,3)=ip6; t(iq1,4)=ip7;
t(iq2,1)=ip1; t(iq2,2)=ip6; t(iq2,3)=ip5; t(iq2,4)=ip7;
t(iq3,1)=ip1; t(iq3,2)=ip2; t(iq3,3)=ip7; t(iq3,4)=ip3;
t(iq4,1)=ip1; t(iq4,2)=ip4; t(iq4,3)=ip3; t(iq4,4)=ip7;
t(iq5,1)=ip1; t(iq5,2)=ip5; t(iq5,3)=ip8; t(iq5,4)=ip7;
t(iq6,1)=ip1; t(iq6,2)=ip4; t(iq6,3)=ip7; t(iq6,4)=ip8;

% boundary nodes
bpx=[ib4 ib2]; bpy=[ib1 ib3]; bpz=[ib5 ib6];
