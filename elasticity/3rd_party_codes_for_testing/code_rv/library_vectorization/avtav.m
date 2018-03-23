function as = avtav (ava,avb)
% ava: ava(1:nx,1,1:nz)
% avb: ama(1:nx,1,1:nz)
% as: as(1:nz,1)

as=squeeze(sum(ava.*avb,1));

return