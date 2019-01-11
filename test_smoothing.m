% load sequence
%orig = load_sequence('~/Remote/avocat/denoising/data/synth/mire-gray/%03d.png', 1, 10);
orig = load_sequence('results/two-iter-smoothing-var-aggregation/orig-%03d.png', 1, 2);
orig = mean(orig,3);
orig0 = orig;
%orig = orig(1:4:end,1:4:end,:,:);

rng('default')
sigma = 20;
nisy = orig + sigma*randn(size(orig));
nisy = load_sequence('../bwd-nlkalman/build/mire-nn-tests/s20/%03d.tif', 1, 3);

%%% % load optical flow
%%% addpath('~/Work/optical_flow/algos/flow-code-matlab/');
%%% bof = readFlowSequence('../../../data/derf/bus/s00/tvl1_%03d_b.flo', 1, 150);

prms1.px = 8;
prms1.wx = 11;
prms1.nx = 60;
prms1.nx_agg = 1;
prms1.nt = 40;
prms1.nt_agg = 1;
prms1.beta_x = 2;
prms1.beta_t = 1;
prms1.r = prms1.px * prms1.px;
prms1.filter_type = 'pos';
prms1.distance = 'L2';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% smoothing with following frame                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filt1  = zeros(size(nisy));
filt2  = zeros(size(nisy));
smoo1  = zeros(size(nisy));

for f = 1:size(nisy,4),

	disp(sprintf('frame % 3d', f))
	if f == 1,
		filt1(:,:,:,f) = bwd_nlkalman_filter(nisy(:,:,:,f), [], [], [], sigma, prms1);
		filt2(:,:,:,f) = filt1(:,:,:,f);

	else
		% filtering of f+1
		filt1(:,:,:,f) = bwd_nlkalman_filter(nisy (:,:,:,f), [], [], ...
		                                     filt2(:,:,:,f-1), sigma, prms1);

		filt2(:,:,:,f) = bwd_nlkalman_filter(nisy(:,:,:,f), [], filt1(:,:,:,f), ...
		                                     filt2(:,:,:,f-1), sigma, prms1);

		% smoothing of f
		smoo1(:,:,:,f-1) = bwd_nlkalman_smoother(filt2(:,:,:,f-1), filt2(:,:,:,f), [],...
		                                         sigma, prms1);

%		imwrite(uint8(smoo1(:,:,1,f-1)),sprintf('smoo1-%03d.png',f-1));
		iio_write(sprintf('smoo1-%03d.tif',f-1), smoo1(:,:,1,f-1));
	end

%	imwrite(uint8(filt2(:,:,1,f)),sprintf('filt2-%03d.png',f));
%	imwrite(uint8(filt1(:,:,1,f)),sprintf('filt1-%03d.png',f));
	iio_write(sprintf('filt1-%03d.tif',f), filt1(:,:,1,f));
	iio_write(sprintf('filt2-%03d.tif',f), filt2(:,:,1,f));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% smoothing with full video                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

smoon  = zeros(size(nisy));
for f = size(nisy,4):-1:1,

	disp(sprintf('frame % 3d', f))
	if f == size(nisy,4),
		smoon(:,:,:,f) = filt2(:,:,:,f);

	else
		smoon(:,:,:,f) = bwd_nlkalman_smoother(filt2(:,:,:,f), smoon(:,:,:,f+1), [], ...
		                                       sigma, prms1);

	end

	imwrite(uint8(smoon(:,:,1,f)),sprintf('smooN-%03d.png',f));
end
