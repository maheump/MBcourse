function [ PofA, PofB ] = EG_EpsilonGreedy( pa, pb, epsilon )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

if numel(pa) ~= numel(pb)
    error('pa and pb must have the same size');
end

N = numel(pa);

Va = pa;
Vb = pb;

[~, stim] = max([Va; Vb], [], 1);
stim = repmat(stim, [numel(epsilon), 1]);

epsilon1 = repmat(epsilon', [1, N]); % for the other choice (exploration)
epsilon2 = repmat(1-epsilon', [1, N]); % for the choice with maximum value

PofA = NaN(numel(epsilon), N);
PofA(stim == 1) = epsilon2(stim == 1);
PofA(stim == 2) = epsilon1(stim == 2);
PofB = NaN(numel(epsilon), N);
PofB(stim == 2) = epsilon2(stim == 2);
PofB(stim == 1) = epsilon1(stim == 1);

end