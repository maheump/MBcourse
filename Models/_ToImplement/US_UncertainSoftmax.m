function [ PofA, PofB ] = US_UncertainSoftmax( pa, pb, phi, beta )
%US_Softmax computes the probability of choosing each action as a function
%of a beta value, a phi value and at each trial (3D matrices).
%   "pa", "pb", "phi" and "beta" must be 1-row arrays.

if numel(pa) ~= numel(pb)
    error('pa and pb must have the same size');
end

N = numel(pa);

Va = US_Vi(pa, phi);
Vb = US_Vi(pb, phi);

Va = reshape(Va, [1, size(Va)]);
Vb = reshape(Vb, [1, size(Vb)]);

Va = repmat(Va, [numel(beta), 1, 1]);
Vb = repmat(Vb, [numel(beta), 1, 1]);

beta = repmat(beta', 1, numel(phi), N);

PofA = exp(beta .* Va) ./ (exp(beta .* Va) + exp(beta .* Vb));
PofB = exp(beta .* Vb) ./ (exp(beta .* Va) + exp(beta .* Vb));

end