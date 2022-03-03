function [PCs,percent_explained]=do_PCA(X);

X = X-mean(X,1);

C_X = cov(X);

[V,D]=eig(C_X);

% Sort the eigenvalues
[eigvals,i]=sort(diag(D),'descend');

percent_explained=100*eigvals/sum(eigvals);

% And re-sort the eigenvectors
V=V(:,i);

% Project to get the "principal components"
PCs=X*V;
end