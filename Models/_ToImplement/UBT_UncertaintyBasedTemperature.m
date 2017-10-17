function [ PofA, PofB ] = UBT_UncertaintyBasedTemperature( pa, pb, beta0 )
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

if numel(pa) ~= numel(pb)
    error('pa and pb must have the same size');
end

N = numel(pa);

Va = repmat(pa, [numel(beta0), 1]);
Vb = repmat(pb, [numel(beta0), 1]);

d = sqrt(mean([pa; pb].^2, 1) - (mean([pa; pb], 1).^2));
d(d < 0.01) = 0.01;
d = repmat(d, [numel(beta0), 1]);

beta0 = repmat(beta0', [1, N]);

betai = beta0 ./ d;

PofA = exp(betai .* Va) ./ (exp(betai .* Va) + exp(betai .* Vb));
PofB = exp(betai .* Vb) ./ (exp(betai .* Va) + exp(betai .* Vb));

end

