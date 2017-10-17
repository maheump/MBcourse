function [ PofA, PofB ] = SM_Softmax( pa, pb, beta )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

if numel(pa) ~= numel(pb)
    error('pa and pb must have the same size');
end

Va = repmat(pa, [numel(beta), 1]);
Vb = repmat(pb, [numel(beta), 1]);

beta = repmat(beta', size(pa));

PofA = exp(beta .* Va) ./ (exp(beta .* Va) + exp(beta .* Vb));
PofB = exp(beta .* Vb) ./ (exp(beta .* Va) + exp(beta .* Vb));

end