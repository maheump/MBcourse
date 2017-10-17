function [ BIC, LSE, gridsearch, bestbeta, bestphi, PofA, PofB ] = US_SimulateUS( choices, rewards, betagrid, phigrid )
%US_SimulateSoftmax finds the best beta and phi parameters of the US
%algorithm given a sequence of binary choices and associated binary
%outcomes (rewards).
%   "choices", "rewards", "betagrid" and "phigrid" must be 1-row arrays.
%   "choices" must be filled with 1s (A) and 2s (B).
%   "rewards" must be filled with 1s (rewards) and 0s (no rewards).

% Get expected reward rates
[pa, pb, N] = UD_GetExpectedReward(choices, rewards);

% 3 dimension matrix:
%	- column: phi,
%   - rows: beta,
%   - slices: trials.
[PofA, PofB] = US_UncertainSoftmax(pa, pb, phigrid, betagrid);

% Find the best parameters values through LSE
Achoices = choices;
Achoices(Achoices == 2) = 0;
diff = PofA - repmat(reshape(Achoices, [1, 1, numel(Achoices)]), [numel(betagrid), numel(phigrid)]);
diff = diff .^ 2;
gridsearch = sum(diff, 3);
LSE = min(gridsearch(:));
[bestbeta, bestphi] = find(gridsearch == LSE);
BIC = log(LSE/N) + (2*log(N));

end