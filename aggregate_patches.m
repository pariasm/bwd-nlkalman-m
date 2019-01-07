% Aggregates a set of n patches of size hxw with ch channels over an image
%
% USAGE: agg_out = aggregate_patches(agg_in, patches, [ph pw], coordinates)
%
%  -> agg_in      : input aggregation image
%  -> patches     : set of patches (h*w*ch x n)
%  -> [ph,pw,p ]  : patch size
%  -> coordinates : (x,y,t) coordinates of top left pixel of each patch (n x 3)
%
%  <- agg_out     : output aggretation image
function agg = aggregate_patches(agg, patches, psz, cc, mode)

	if (nargin < 5),
		mode = 'avg';
	end

	ph = psz(1); pw = psz(2); pt = psz(3);

	% number of patches
	n = size(cc,1);

	patches = reshape(patches(:,1:n), [ph pw size(agg,3) pt n]);

	if strcmp(mode, 'min'),

		for i = 1:n,

			agg(cc(i,2):cc(i,2)+ph-1,...
				 cc(i,1):cc(i,1)+pw-1,:,...
				 cc(i,3):cc(i,3)+pt-1) = min(patches(:,:,:,i),...
				                           agg(cc(i,2):cc(i,2)+ph-1,...
				                               cc(i,1):cc(i,1)+pw-1,:,...
				                               cc(i,3):cc(i,3)+pt-1));

		end

	else,

		for i = 1:n,

			agg(cc(i,2):cc(i,2)+ph-1,...
				 cc(i,1):cc(i,1)+pw-1,:,...
				 cc(i,3):cc(i,3)+pt-1) = patches(:,:,:,i) ...
				                         + agg(cc(i,2):cc(i,2)+ph-1,...
				                               cc(i,1):cc(i,1)+pw-1,:,...
				                               cc(i,3):cc(i,3)+pt-1);

		end

	end
end

