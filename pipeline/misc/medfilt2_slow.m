function f = medfilt2_slow(f, wsize)

if mod(wsize, 2) ~= 1
	error 'Filter window size must be an odd number.';
end

wr = (wsize - 1) / 2;   % Filter window radius

padded = zeros(size(f, 1) + 2*wr(1), size(f, 2) + 2*wr(2));
padded(wr(1)+1:wr(1)+size(f, 1), wr(2)+1:wr(2)+size(f, 2)) = f;

out_padded = zeros(size(padded));

for r = wr(1)+1:wr(1)+size(f, 1)
	for c = wr(2)+1:wr(2)+size(f, 2)
		window = padded(r-wr(1):r+wr(1), c-wr(2):c+wr(2));
		out_padded(r, c) = median(window(:));
	end
end

f = out_padded(wr(1)+1:wr(1)+size(f, 1), wr(2)+1:wr(2)+size(f, 2));

