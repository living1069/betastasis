function x = unique_preserve_order(x)

[x, idx] = unique(x, 'first');
[~, orig_order] = sort(idx);
x = x(orig_order);

