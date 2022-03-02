function [DLE,who]=get_DLE(X,gridLF,grid_idx)
[where,who]=max(sum(X.^2,2));
tmp=gridLF.pos(gridLF.inside,:);
DLE=pdist2(tmp(who,:),tmp(grid_idx,:));
end
