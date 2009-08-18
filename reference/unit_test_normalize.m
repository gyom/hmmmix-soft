fprintf('Doing a quick unit term to make sure that we have\na function normalize.m that behaves as expected.\n');

A = rand(10,10);

E = normalize(A,1);
for i=1:10
    assert(all(abs(sum(E,1) - 1) < 1e-8));
end

E = normalize(A,2);
for i=1:10
    assert(all(abs(sum(E,2) - 1) < 1e-8));
end

E = normalize(A);
assert(abs(sum(E(:)) - 1) < 1e-8);

fprintf('\tnormalize.m behaves as expected.\n');