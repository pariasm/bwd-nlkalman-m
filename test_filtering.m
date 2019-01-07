% load sequence
orig = load_sequence('~/Remote/avocat/denoising/data/synth/mire-gray/%03d.png', 1, 10);
orig = mean(orig,3);
orig0 = orig;
orig = orig(1:4:end,1:4:end,:,:);


%%% % load optical flow
%%% addpath('~/Work/optical_flow/algos/flow-code-matlab/');
%%% bof = readFlowSequence('../../../data/derf/bus/s00/tvl1_%03d_b.flo', 1, 150);

prms1.px = 8;
prms1.wx = 21;
prms1.nx = 60;
prms1.nx_agg = 60;
prms1.nt = 40;
prms1.nt_agg = 2;
prms1.beta_x = 2;
prms1.beta_t = 2;
prms1.r = prms1.px * prms1.px;
prms1.filter_type = 'pos';
prms1.distance = 'L2';

rng('default')
sigma = 20;
nisy = orig + sigma*randn(size(orig));

deno1 = zeros(size(nisy));
deno2 = zeros(size(nisy));

for f = 1:size(nisy,4),

	disp(sprintf('frame % 3d', f))
	if f == 1,
		deno2(:,:,:,f) = bwd_nlkalman_filter(nisy(:,:,:,f  ), [], [], [], sigma, prms1);
	else
		deno1(:,:,:,f) = bwd_nlkalman_filter(nisy(:,:,:,f  ), [], [], ...
		                                     deno2(:,:,:,f-1), sigma, prms1);

		deno2(:,:,:,f) = bwd_nlkalman_filter(nisy (:,:,:,f  ), [], ...
		                                     deno1(:,:,:,f  ), deno2(:,:,:,f-1), ...
		                                     sigma, prms1);
	end
	
	imwrite(uint8(deno1(:,:,1,f)),sprintf('deno1-%03d.png',f));
	imwrite(uint8(deno2(:,:,1,f)),sprintf('deno2-%03d.png',f));
end


