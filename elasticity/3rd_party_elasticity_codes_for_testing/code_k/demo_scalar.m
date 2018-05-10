%----------------------------------------
% Scalar problem 2D/3D
% -\Delta u=f in \Omega=(0,1)^d
%         u=0 on boundary
% 2D: f=-2x(x-1)-2y(y-1);
%     u_{exact}= xy(x-1)(y-1)
% 3D: f=-2xy(x-1)(y-1)-2yz(y-1)(z-1)-2xz(x-1)(z-1);
%     u_{exact}= xyz(x-1)(y-1)(z-1)
%----------------------------------------
%
d=3;         % dimension
nx=9;        % number of points in each dimension
% mesh generation
if (d==2)
    [p,t,pbx,pby]=kpde2dumsh(0,1,0,1,nx,nx);
    ibcx=union(pbx(:,1),pbx(:,2)); ibcy= union(pby(:,1),pby(:,2));
    ibcd=union(ibcx,ibcy);   
    clear ibcx ibcy
elseif (d==3)
    [p,t,pbx,pby,pbz]=kpde3dumsh(0,1,0,1,0,1,nx,nx,nx);
    ibcx=union(pbx(:,1),pbx(:,2));  ibcy=union(pby(:,1),pby(:,2));
    ibcz=union(pbz(:,1),pbz(:,2)); ibcd=union(ibcx,union(ibcy,ibcz));  
    clear ibcx ibcy ibcz
else
    error('The value of d is not valid')
end

np=size(p,1); nt=size(t,1);
fprintf('Dimension             :%3d\n',d)
fprintf('Number of vertices    :%5d\n',np)
fprintf('Number of elements    :%5d\n',nt)
 
% Right-hand side / exact solution
if (d==2)
   x=p(:,1); y=p(:,2);  ue=x.*y.*(x-1).*(y-1);
   fh=-2*x.*(x-1)-2*y.*(y-1);
   clear x y
else 
    x=p(:,1); y=p(:,2);  z=p(:,3); ue=x.*y.*z.*(x-1).*(y-1).*(z-1);
   fh=-2*x.*y.*(x-1).*(y-1)-2*y.*z.*(y-1).*(z-1)-2*x.*z.*(x-1).*(z-1);
   clear x y z
end

% Assembly of matrices and rhs
if (d==2)
    K=kpde2dstf(p,t,1);
    M=kpde2dmss(p,t,1);
    f=kpde2drhs(p,t,fh);
else
    K=kpde3dstf(p,t,1);
    M=kpde3dmss(p,t,1);
    f=kpde3drhs(p,t,fh);
end

% Boundary conditions (penalization)
K(ibcd,ibcd)=K(ibcd,ibcd)+10^15*speye(length(ibcd));
f(ibcd)=0;
 
% Solution (Gaussian elimination)
u=K\f;

% Visualisation
if (d==2)
   trisurf(t,p(:,1),p(:,2),u,'facecolor','interp')
   grid off, colorbar('horiz')
else
    pt=(p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:)+p(t(:,4),:))/4;
   it=find(pt(:,2)>=1/2);
   kpde3dshow(p,t,u,it)
end

% Approximate L^2-error
du=u-ue; err=sqrt(du'*M*du);
fprintf('Approximate L^2-error : %12.8e \n',err)
