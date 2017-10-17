function [ BIC, LSE, gridsearch, bestbeta, PofA, PofB ] = UBT_SimulateUBT( choices, rewards, beta0grid )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

% Get expected reward rates
[pa, pb, N] = UD_GetExpectedReward(choices, rewards);

% 2 dimension matrix:
%	- trials: beta0,
%   - rows: beta.
[PofA, PofB] = UBT_UncertaintyBasedTemperature(pa, pb, beta0grid);

% Find the best parameters values through LSE
Achoices = choices;
Achoices(Achoices == 2) = 0;
diff = PofA - repmat(Achoices, [numel(beta0grid), 1]);
diff = diff .^ 2;
gridsearch = sum(diff, 2);
LSE = min(gridsearch(:));
bestbeta = find(gridsearch == LSE);
BIC = log(LSE/N) + (1*log(N));

end

