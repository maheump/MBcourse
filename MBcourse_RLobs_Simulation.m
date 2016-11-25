function out = MBcourse_RLobs_Simulation(alpha, linkSpec, fullOutcome, nSimus, subchoices)

% Complete output
if nargin < 4, nSimus = 1; end
if ~iscell(linkSpec), linkSpec = {linkSpec}; end

% Get the task information
[nTrials, nStims] = size(fullOutcome);

% Prepare outputs
choice        = NaN(nTrials, nSimus);
choiceOutcome = NaN(nTrials, nSimus);
values        = NaN(nTrials, nStims, nSimus);
choiceProba   = NaN(nTrials, nStims, nSimus);

% For each agent to simulate
for s = 1:nSimus
    
    % Define initial value of each stimulus
    v0 = ones(1,2) ./ nStims;
    v  = v0;
    
    % For each trial
    for t = 1:nTrials
        
        % Compute likelihood of each choice option
        p = MBcourse_RLobs_ObsFuns(v, linkSpec);
        
        % Do a weighted coinflip to make a choice:
        %   - choose stim 1 if a random number is in the [0, p(1)[ interval,
        %   - choose stim 2 if a random number is in the [p(1), 1] interval,
        if nargin < 5
            if rand(1) < p(1), c = 1;
            else,              c = 2; end
            
        % Or use choices already made by the subject (used for the fitting
        % procedure
        else, c = subchoices(t);
        end
        choice(t,s) = c;
        
        % Select the outcome (i.e. the reward) for this choice
        r = fullOutcome(t,c);
        
        % Store stimuli values and choice probability
        values(t,:,s) = v; % store value of V
        choiceProba(t,:,s) = p;
        
        % Update values
        pe  = r - v(c); % compute prediction error
        v(c)= v(c) + (alpha * pe); % update value
    end
    
    % Store the feedback the subject received
    c1 = choice(:,s) == 1;
    c2 = choice(:,s) == 2;
    choiceOutcome(c1,s) = fullOutcome(c1,1);
    choiceOutcome(c2,s) = fullOutcome(c2,2);
end

% Export the information
out               = [];
out.alpha         = alpha;
out.linkSpec      = linkSpec;
out.choice        = choice;
out.choiceOutcome = choiceOutcome;
out.values        = values;
out.choiceProba   = choiceProba;

end