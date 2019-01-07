% cc = [x, y] = [col, row]
function patches = extract_patches(cc, psz, im)

h = size(im,1);
w = size(im,2);
c = size(im,3);
n = size(cc,1);
patches = nan(psz*psz*c,n);

cc = round(cc);
in = find(cc(:,1) >= 1 & cc(:,1) <= w - psz + 1 & ...
          cc(:,2) >= 1 & cc(:,2) <= h - psz + 1);

for i = in',
	x = cc(i,1);
	y = cc(i,2);
	patches(:,i) = reshape(im(y:y+psz-1,x:x+psz-1,:), [psz*psz*c 1]);
end



