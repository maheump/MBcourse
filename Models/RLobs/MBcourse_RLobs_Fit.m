function [ MSE, LLH ] = MBcourse_RLobs_Fit( pChoice, subChoice )

% Get trials' list associated to each choice
c1 = subChoice == 1; % A
c2 = subChoice == 2; % B

% Compute the mean squared error
MSE = mean([(pChoice(c1,1) - 1) .^2; ...
            (pChoice(c2,2) - 1) .^2]);

% Compute the log likelihood with the p(choice|model)
LLH = sum(log( pChoice(c1, 1))) ...
    + sum(log( pChoice(c2, 2)));

end