function [p,t,bpx,bpy]=kpde2dumsh(ax,bx,ay,by,nx,ny)
%KPDE2DUMSH Uniform mesh generation of (ax,bx)*(ay,by)
%--------------------------------------------------------------------
%      [p,t,bpx,bpy]=kpde2dumsh(ax,bx,ay,by,nx,ny)
% 
% Inputs:
%   (ax,bx) : x-interval
%   (ay,by) : y-interval
%     nx,ny : number of points in x and y
%
% Ouputs:
%        p  : nodes coordinates, array np*2
%        t  : triangle vertices, array nt*3
%  bpx,bpy  : boundary nodes, arrays nx*2 and ny*2
%             bpx(:,1), bpx(:,2), bundary nodes at y=ay and y=by 
%             bpy(:,1), bpy(:,2), boundary nodes at x=ax and x=bx
%--------------------------------------------------------------------
np=nx*ny; nq=(nx-1)*(ny-1); nt=2*nq;               

% vertce coordinates
hx=(bx-ax)/(nx-1); hy=(by-ay)/(ny-1);
[x,y]=ndgrid(ax:hx:bx,ay:hy:by);
p=[x(:),y(:)];

% quadrangles
ip=[1:nx*ny]';
ib1=[1:nx]'; ib2=nx*[1:ny]';
ib3=[nx*(ny-1)+1:nx*ny]'; ib4=[1:nx:nx*ny]';
ib23=union(ib2,ib3); ib34=union(ib3,ib4);
ib14=union(ib1,ib4); ib12=union(ib1,ib2);
iq1=setdiff(ip,ib23); iq2=setdiff(ip,ib34);
iq3=setdiff(ip,ib14); iq4=setdiff(ip,ib12);

% triangulation
t=zeros(nt,3);
t(1:nq,1)=iq1;      t(1:nq,2)=iq2;      t(1:nq,3)=iq3;
t(nq+1:2*nq,1)=iq3; t(nq+1:2*nq,2)=iq4; t(nq+1:2*nq,3)=iq1;

% boundary nodes
bpx=[ib1 ib3];
bpy=[ib4 ib2];