function [ BIC, LSE, gridsearch, bestepsilon, PofA, PofB ] = EG_SimulateEG( choices, rewards, epsilongrid )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% Get expected reward rates
[pa, pb, N] = UD_GetExpectedReward(choices, rewards);

% 3 dimension matrix:
%	- column: trials,
%   - rows: epsilon.
[PofA, PofB] = EG_EpsilonGreedy(pa, pb, epsilongrid);

% Find the best parameters values through LSE
Achoices = choices;
Achoices(Achoices == 2) = 0;
diff = PofA - repmat(Achoices, [numel(epsilongrid), 1]);
diff = diff .^ 2;
gridsearch = sum(diff, 2)';
LSE = min(gridsearch(:));
bestepsilon = find(gridsearch == LSE);
BIC = log(LSE/N) + (1*log(N));

end