function p = fexact_wizard(a, b)

a = a(:); b = b(:);

fprintf('\t+A\t-A\tTOT\n');
fprintf('+B\t%d\t%d\t%d\n', sum(a & b), sum(~a & b), sum(b));
fprintf('-B\t%d\t%d\t%d\n', sum(a & ~b), sum(~a & ~b), sum(~b));
fprintf('TOT\t%d\t%d\t%d\n', sum(a), sum(~a), length(a));

p = fexact(a, b);

