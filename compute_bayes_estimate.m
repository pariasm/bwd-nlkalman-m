function [deno, var] = ...
                 compute_bayes_estimate(nisy, bsic, sigma, beta, rank, filter_type, mu, U)

if nargin < 5,
	filter_type = 'pos';
end

if nargin < 6,
	mu = [];
	U = [];
end

sigma_bsic = 0;
if isempty(bsic),
	bsic = nisy;
	sigma_bsic = sigma;
end

d  = size(bsic,1);
nb = size(bsic,2);
nn = size(nisy,2);

% center
if isempty(mu), mean_bsic = mean(bsic, 2);
else            mean_bsic = mu;
end
mean_nisy = mean_bsic;

nisy = nisy - mean_nisy*ones(1,nn);
bsic = bsic - mean_bsic*ones(1,nb);

if isempty(U),
	% redefine rank as...
	if (rank >= 0), rank = min(rank, nb-1);
	else            rank = min(d   , nb-1);
	end

	% eigendecomposition
	if rank > 0 && rank < d-1,
		[U,S,V] = svds(bsic,rank);
	else
		[U,S,V] = svd(bsic);
	end

	% turn singular values of data matrix into cov. eigenvectors
	S = diag(S).^2/nb;

else
	% project data vectors onto basis
	S = mean((U'*bsic).^2,2);
end

% sort basis according to variance
if isempty(U) && rank > 0 && rank <= d,
	tmp = sortrows([ S U'],-1);
	S = tmp(:,1);
	U = tmp(:,2:end)';
end

% wiener filter coefficients
if strcmp(filter_type, 'pos'),

	S = min(max( 0 , S - sigma_bsic*sigma_bsic),Inf);
	W = diag(1./(1 + beta * sigma*sigma ./ S));

elseif strcmp(filter_type, 'pca'),

	W = diag(ones(size(S)));

elseif strcmp(filter_type, 'neg'),

	S = min(max(-Inf, S - sigma_bsic*sigma_bsic),Inf);
	W = diag(1./(1 + beta * sigma*sigma ./ S));

elseif strcmp(filter_type, 'neg-inv'),

	S = min(max(-Inf, S - sigma_bsic*sigma_bsic),Inf);
	W = abs(diag(1./(1 + beta * sigma*sigma ./ S)));

end

% denoise group
comp = U'*nisy;
var = S.*diag(W);
deno = (U*W)*comp + mean_nisy*ones(1,nn);

