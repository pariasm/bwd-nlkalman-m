% Backward NL-Kalman filter. 
%
% USAGE: deno = bwd_nlkalman_filter(nisy1, nisy0, filt1, filt0, sigma, prams)
%
%  -> nisy1    : noisy frame t
%  -> nisy0    : noisy frame t-1 (optional for parameter estimation)
%  -> filt1    : pilot estimate of frame t (optional for parameter estimation)
%  -> filt0    : denoised frame t-1 (leave empty for the first frame)
%  -> sigma    : noise std. dev.
%  -> prams    : struc with prams (wx, px, np, r, lambda)
%
%  <- deno     : denoised image
function [recon, recon_x, recon_t, aggw_x, aggw_t] = ...
                          bwd_nlkalman_filter(nisy1, nisy0, deno1, deno0, sigma, prams)

% image size and channels
w = size(nisy1,2);
h = size(nisy1,1);
chnls = size(nisy1,3);

% aggregation
aggw   = zeros(size(nisy1)); % aggregation weights
aggw_t = zeros(size(nisy1)); % aggregation weights
aggw_x = zeros(size(nisy1)); % aggregation weights
aggu_t = zeros(size(nisy1)); % aggregated image
aggu_x = zeros(size(nisy1)); % aggregated image

recon   = zeros(size(nisy1)); % denoised image (recon = aggu ./ aggw);
recon_t = zeros(size(nisy1)); % denoised image (recon = aggu ./ aggw);
recon_x = zeros(size(nisy1)); % denoised image (recon = aggu ./ aggw);

% patch dimensionality
pdim = prams.px * prams.px * chnls;


% vars0  = zeros(h,w,pdim);
% vars01 = zeros(h,w,pdim);


% step sizes
psz = prams.px;
wsz = prams.wx;
stepx = max(1, floor(psz/2));
stepy = max(1, floor(psz/2));
%stepx = 1;
%stepy = 1;
jj = 0; % iteration counter (only for display)

% dct basis
U = dct_basis(psz,psz);

for pay = [1 : stepy : h - psz + 1, h - psz + 1],
for pax = [1 : stepx : w - psz + 1, w - psz + 1], jj = jj + 1;

	% process only if a pixel inside the patch has low aggregation weight
	aggw_min = min(reshape(aggw(pay:pay + psz - 1, pax:pax + psz - 1,1),[pdim,1]));
	if aggw_min > 1, continue; end

	% -------------------------------------------------- compute patch group
	% patches in search region
	srx = [max(1,pax - wsz): min(w ,pax + wsz + psz - 1)];
	sry = [max(1,pay - wsz): min(h ,pay + wsz + psz - 1)];

	if isempty(deno1),
		srch_region = nisy1(sry, srx, :);
		refe_patch  = nisy1(pay:pay + psz - 1, pax:pax + psz - 1,:);
	else
		srch_region = deno1(sry, srx, :);
		refe_patch  = deno1(pay:pay + psz - 1, pax:pax + psz - 1,:);
	end
	srch_patches = im2col_ch(srch_region, [psz psz]);

	no_time = sum(isnan(refe_patch(:))) | isempty(deno0);

	if strcmp(prams.distance, 'L2'),
		patches.distas = L2_distance(refe_patch(:), srch_patches);
	else
		patches.distas = L1_distance(refe_patch(:), srch_patches);
	end

	% number of patches
	if no_time, npatches = min(prams.nx, length(patches.distas));
	else        npatches = min(prams.nt, length(patches.distas));
	end

	% sort distances
	[patches.distas, idx] = sort(patches.distas);

	% coordinates of the np nearest neighbors to the ref patch
	idx = idx(1:npatches)';

	srch_h = size(srch_region,1) - psz + 1;
	patches.coords = [max(1,pax - prams.wx) + floor((idx-1)/srch_h),...
	                  max(1,pay - prams.wx) +   mod( idx-1 ,srch_h), ones(npatches,1)];

	% ---------------------------------------------- extract similar patches
	Y1 = double(srch_patches(:,idx));

	% ----------------------------------------------- compute bayes estimate
	if no_time, npatches_agg = min(prams.nx_agg, npatches);
	else        npatches_agg = min(prams.nt_agg, npatches);
	end
	gg = ones(pdim, npatches_agg); % aggregation weights are 1

	% attemp temporal denoising
	if no_time == 0,

		% learn state covariance
		srch_region = deno0(sry, srx, :);
		srch_patches_deno0 = im2col_ch(srch_region, [psz psz]);
		X0 = U' * srch_patches_deno0(:,idx);

		idx_nans = find(sum(isnan(X0)));
		X0(:,idx_nans) = [];
		Y1(:,idx_nans) = [];
		patches.coords(idx_nans,:) = [];

		if ~isempty(X0),
			var0 = var(X0,1,2);

			% learn transition variances
			use_clean = isempty(nisy0);
			Y1 = U' * Y1;
			if use_clean,
				if ~isempty(deno1),
					srch_region = deno1(sry, srx, :);
					srch_patches_deno1 = im2col_ch(srch_region, [psz psz]);
					X1 = U' * srch_patches_deno1(:,idx);
					X1(:,idx_nans) = [];
					var01 = max(0, mean((X1 - X0).^2, 2));
				else
					var01 = max(0, mean((Y1 - X0).^2, 2) - sigma^2);
				end
			else
				srch_region = nisy0(sry, srx, :);
				srch_patches_nisy0 = im2col_ch(srch_region, [psz psz]);
				Y0 = U' * srch_patches_nisy0(:,idx);
				Y0(:,idx_nans) = [];
				var01 = max(0, mean((Y1 - Y0).^2, 2) - 2*sigma^2);
			end

			% Kalman filter
			X0 = mean(X0(:,1:npatches_agg),2)*ones(1,npatches_agg);
			if isempty(deno1),
				Y1 = Y1(:,1:npatches_agg);
			else
				srch_region = nisy1(sry, srx, :);
				srch_patches_nisy1 = im2col_ch(srch_region, [psz psz]);
				Y1 = srch_patches_nisy1(:,idx);
				Y1(:,idx_nans) = [];
				Y1 = U' * Y1(:,1:npatches_agg);
			end

			var1 = var0 + var01;
			K = var1./(var1 + prams.beta_t * sigma^2);   % Kalman gains
			X1 = U * ( X0 + diag(K)*(Y1 - X0) );

			% variances
			var1 = (1-K).^2.*var1 + K.^2*sigma^2;

			% aggregate
			cc = patches.coords(1:npatches_agg,:);
			ivar1 = 1/max(sum(var1), 1e-4);
			aggu_t = aggregate_patches(aggu_t, ivar1*gg.*X1, [psz psz 1], cc);
			aggw_t = aggregate_patches(aggw_t, ivar1*gg    , [psz psz 1], cc);

%			vars0 (pay, pax, :) = var0;
%			vars01(pay, pax, :) = var01;


		else
			no_time = 1;

		end

	end

	if no_time
		% spatial denoising
		r = prams.r;
		beta = prams.beta_x;
		type = prams.filter_type;
		[X1, var1] = compute_bayes_estimate(Y1, [], sigma, beta, r, type, [], U);

		cc = patches.coords(1:prams.nx_agg,:);
		ivar1 = 1/max(sum(var1), 1e-6);
		aggu_x = aggregate_patches(aggu_x, ivar1*gg.*X1, [psz psz 1], cc);
		aggw_x = aggregate_patches(aggw_x, ivar1*gg    , [psz psz 1], cc);

	end

	% ------------------------------------------- aggregate patches on image
	aggu = aggu_t + aggu_x;
	aggw = aggw_t + aggw_x;


	% draw
	if (mod(jj,100) == 1)
		nz = find(aggw ~= 0);
		recon(nz) = aggu(nz) ./ aggw(nz);
		imagesc(max(min([recon],255),0)/255,[0 1]);
		axis equal, axis off, colormap gray, drawnow, pause(.01)
	end

end
end

% compute denoised image
nonzero = find(aggw ~= 0);
%recon(nonzero) = min(255, max(0, aggu(nonzero) ./ aggw(nonzero)));
recon(nonzero) = aggu(nonzero) ./ aggw(nonzero);

nonzero = find(aggw_t ~= 0);
recon_t(nonzero) = aggu_t(nonzero) ./ aggw_t(nonzero);

nonzero = find(aggw_x ~= 0);
recon_x(nonzero) = aggu_x(nonzero) ./ aggw_x(nonzero);
