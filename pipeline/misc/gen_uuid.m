function uuid = gen_uuid()
persistent rand_seeded_;
if isempty(rand_seeded_)
	RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));
	rand_seeded_ = true;
end

uuid = '';
for k = 1:4
	uuid = [uuid dec2hex(randi(2^32 - 1, 1), 8)];
end
