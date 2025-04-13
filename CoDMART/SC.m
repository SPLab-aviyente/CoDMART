  function [U] = SC(W, k)

% calculate degree matrix
degs = sum(W, 2);
degs=degs.^(-0.5);
degs(degs==inf)=0;
degs=diag(degs);
% D    = sparse(1:size(W, 1), 1:size(W, 2), degs);

% compute unnormalized Laplacian
% L = D - W;
L = W;


% compute normalized Adjacency

% avoid dividing by zero
%         degs(degs == 0) = eps;
        % calculate D^(-1/2)
%         D = spdiags(1./(degs.^0.5), 0, size(D, 1), size(D, 2));
%         D(D==inf)=0;
        % calculate normalized Adjacency
        L = degs * L * degs;

        L=(L+L')/2;

% compute the eigenvectors corresponding to the k largest
% eigenvalues
% diff   = eps;
% [U, VV] = eigs(L,k,'largestabs');
       [UU, VV] = eig(L);
       v=diag(VV);
       [~, is] = sort(v,'descend');
	   EVec = UU(:,is);
 
       U=normc(EVec(:,1:k));


end