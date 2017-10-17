function [ BIC, LSE, gridsearch, bestalpha, bestbeta, PofA, PofB ] = RL_SimulateStandardRL( choices, rewards, alphagrid, betagrid )
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here

N = numel(choices);

% Run the Rescorla-Wagner algorithm
[Va, Vb] = RL_GetOptionsValues(choices, rewards, alphagrid);

% 2 dimension matrix:
%	- column: trials,
%   - rows: beta,
%   - slices: alpha.
Va = reshape(Va', [1, N, numel(alphagrid)]);
Vb = reshape(Vb', [1, N, numel(alphagrid)]);

% Apply a classical softmax
[PofA, PofB] = SM_Softmax(Va, Vb, betagrid);

% Find the best parameters values through LSE
Achoices = choices;
Achoices(Achoices == 2) = 0;
diff = PofA - repmat(Achoices, [numel(betagrid), 1, numel(alphagrid)]);
diff = diff .^ 2;
gridsearch = squeeze(sum(diff, 2));
LSE = min(gridsearch(:));
[bestbeta, bestalpha] = find(gridsearch == LSE);
BIC = log(LSE/N) + (2*log(N));

end