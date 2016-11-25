function out = MBcourse_GenerateTaskDesign( nTrialsPerBlocks, pReward )

% Get the total number of trials
nTrials = sum(nTrialsPerBlocks);
nBlocks = numel(nTrialsPerBlocks);

% Get reward probability for one of the stimuli in each block
probTrial = repmat([pReward, 1-pReward], 1, nBlocks/2);

% Prepare outputs
feedbackprob = NaN(nTrials, 1);
feedback     = NaN(nTrials, 2);

% Generate the feedback sequence
t = 1;
for b = 1:nBlocks
    for x = 1:nTrialsPerBlocks(b)
        feedbackprob(t) = probTrial(b);
        feedback(t,1) = double(rand(1) <= feedbackprob(t));
        feedback(t,2) = 1 - feedback(t,1);
        t = t+1;
    end
end

% Export all that information 
out = struct('nTrialsPerBlocks', nTrialsPerBlocks, 'pReward', pReward, ...
             'nTrials', nTrials, 'nBlocks', nBlocks, 'probTrial', probTrial, ...
             'feedbackprob', feedbackprob, 'feedback', feedback);

end