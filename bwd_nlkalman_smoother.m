% Backward NL-Kalman smoother.
%
% USAGE: deno = bwd_nlkalman_smoother(filt1, smoo0, smoo1, sigma, prams)
%
%  -> filtered1 : Kalman filtered image at frame t
%  -> smoothed0 : Kalman smoothed image at frame t+1 (leave empty for the last frame)
%  -> smoothed1 : Kalman smoothed image at frame t (optional for parameter estimation)
%  -> sigma     : noise std. dev.
%  -> prams     : struc with parameters
%
%  <- deno      : smoothed image image at frame t
function [recon, recon_x, recon_t, aggw_x, aggw_t] = ...
                 bwd_nlkalman_smoother(filtered1, smoothed0, smoothed1, sigma, prams)

% image size and channels
w = size(filtered1,2);
h = size(filtered1,1);
chnls = size(filtered1,3);

% aggregation
aggw   = zeros(size(filtered1)); % aggregation weights
aggw_t = zeros(size(filtered1)); % aggregation weights
aggw_x = zeros(size(filtered1)); % aggregation weights
aggu_t = zeros(size(filtered1)); % aggregated image
aggu_x = zeros(size(filtered1)); % aggregated image

recon   = zeros(size(filtered1)); % denoised image (recon = aggu ./ aggw);
recon_t = zeros(size(filtered1)); % denoised image (recon = aggu ./ aggw);
recon_x = zeros(size(filtered1)); % denoised image (recon = aggu ./ aggw);

% patch dimensionality
pdim = prams.px * prams.px * chnls;


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

	% reference patch
	if isempty(smoothed1),
		refe_patch = filtered1(pay:pay + psz - 1, pax:pax + psz - 1,:);
	else
		refe_patch = smoothed1(pay:pay + psz - 1, pax:pax + psz - 1,:);
	end

	no_time = sum(isnan(refe_patch(:))) | isempty(smoothed0);

	% attemp temporal denoising
	if no_time == 0,

		% -------------------------------------------------- compute patch group
		% patches in search region
		srx = [max(1,pax - wsz): min(w ,pax + wsz + psz - 1)];
		sry = [max(1,pay - wsz): min(h ,pay + wsz + psz - 1)];

		if isempty(smoothed1), srch_region = filtered1(sry, srx, :);
		else                   srch_region = smoothed1(sry, srx, :);
		end
		srch_patches = im2col_ch(srch_region, [psz psz]);

		if strcmp(prams.distance, 'L2'),
			patches.distas = L2_distance(refe_patch(:), srch_patches);
		else
			patches.distas = L1_distance(refe_patch(:), srch_patches);
		end

		% number of patches
		npatches = min(prams.nt, length(patches.distas));

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
		npatches_agg = min(prams.nt_agg, npatches);
		gg = ones(pdim, npatches_agg); % aggregation weights are 1

		% learn state covariance
		srch_region = smoothed0(sry, srx, :);
		srch_patches_smoothed0 = im2col_ch(srch_region, [psz psz]);
		X0 = U' * srch_patches_smoothed0(:,idx);

		idx_nans = find(sum(isnan(X0)));
		X0(:,idx_nans) = [];
		patches.coords(idx_nans,:) = [];
		var0 = var(X0,1,2);

		if ~isempty(X0),

			% learn transition variances
			Y1 = U' * Y1;
			var1 = var(Y1,1,2);
			if ~isempty(smoothed1),
				srch_region = smoothed1(sry, srx, :);
				srch_patches_smoothed1 = im2col_ch(srch_region, [psz psz]);
				X1 = U' * srch_patches_smoothed1(:,idx);
				X1(:,idx_nans) = [];
				var01 = mean((X1 - X0).^2, 2);
			else
				var01 = mean((Y1 - X0).^2, 2);
			end

			% Kalman filter
			X0 = mean(X0(:,1:npatches_agg),2)*ones(1,npatches_agg);
			if isempty(smoothed1),
				Y1 = Y1(:,1:npatches_agg);
			else
				srch_region = filtered1(sry, srx, :);
				srch_patches_filtered1 = im2col_ch(srch_region, [psz psz]);
				Y1 = srch_patches_filtered1(:,idx);
				Y1(:,idx_nans) = [];
				Y1 = U' * Y1(:,1:npatches_agg);
			end

			K = var1./(var1 + var01);   % Kalman gains
			X1 = U * ( Y1 + diag(K)*(X0 - Y1) );

			% variances
			var1 = var1 + K.^2.*max(var0 - var1 - var01, 0);

			% aggregate
			cc = patches.coords(1:npatches_agg,:);
			ivar1 = 1/max(sum(var1), 1e-4);
			aggu_t = aggregate_patches(aggu_t, ivar1*gg.*X1, [psz psz 1], cc);
			aggw_t = aggregate_patches(aggw_t, ivar1*gg    , [psz psz 1], cc);

		else
			no_time = 1;

		end

	end

	if no_time,
		% we just copy the previous filtered patch with very small agg weights
		X1 = refe_patch;

		ivar1 = 1e-6;
		aggu_x = aggregate_patches(aggu_x, ivar1*X1, [psz psz 1], [pax, pay, 1]);
		aggw_x = aggregate_patches(aggw_x, ivar1   , [psz psz 1], [pax, pay, 1]);

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
