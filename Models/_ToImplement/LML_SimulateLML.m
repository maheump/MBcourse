function [ BIC, LSE, PofA, PofB ] = LML_SimulateLML( choices, rewards )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Get expected reward rates
[pa, pb, N] = UD_GetExpectedReward(choices, rewards);

% In the local matching law, the probability to choose an action i
% (amongst two rewarding location) is given by:
[PofA, PofB] = LML_LocalMatchingLaw(pa, pb);

% Compute the error (there is no free parameters)
Achoices = choices;
Achoices(Achoices == 2) = 0;
diff = PofA - Achoices;
diff = diff .^ 2;
LSE = sum(diff);
BIC = log(LSE/N) + (0*log(N));

end