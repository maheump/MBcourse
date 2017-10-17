function [ BIC, LSE, gridsearch, bestbeta, PofA, PofB ] = SM_SimulateSoftmax( choices, rewards, betagrid )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% Get expected reward rates
[pa, pb, N] = UD_GetExpectedReward(choices, rewards);

% 2 dimension matrix:
%	- column: trials,
%   - rows: beta.
[PofA, PofB] = SM_Softmax(pa, pb, betagrid);

% Find the best parameters values through LSE
Achoices = choices;
Achoices(Achoices == 2) = 0;
diff = PofA - repmat(Achoices, [numel(betagrid), 1]);
diff = diff .^ 2;
gridsearch = sum(diff, 2);
LSE = min(gridsearch(:));
bestbeta = find(gridsearch == LSE);
BIC = log(LSE/N) + (1*log(N));

end