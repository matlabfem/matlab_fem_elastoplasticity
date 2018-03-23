function avb = svamt (svx,ama)
% ama: ama(1:ny,1:nx,1:nz)
% svx: svx(1,1:nx)
% avb: avb(1,1:ny,1:nz)

[ny,nx,nz] = size(ama);

avx = svx;
avx = avx(ones(ny,1),:,ones(nz,1));

avb = ama .* avx;
avb = sum(avb,2);
avb = reshape(avb,1,ny,nz);

return