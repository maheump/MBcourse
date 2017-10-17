function p = MBcourse_RLobs_ObsFuns( values, linkspec )

% Make sure that "linspec" is a cell array
if ~iscell(linkspec), linkspec = {linkspec}; end

% Get the link function
fun = linkspec{1};

% Local matching rule
% ===================
if strcmpi(fun, 'LMR')
    
    % Compute the sum of the values
    sev = sum(values);
    
    % Divide the value of each options by this sum such that we have a
    % probability of choosing each action
    p = values ./ sev;
    
% Epsilon-greedy choice rule
% ==========================
elseif strcmpi(fun, 'Greedy')
	
    % Free parameter of the epsilon-greedy choice rule
    epsilon = linkspec{2};
    
    % The best option is the one with the maximum value
    [~, bestoption] = max(values);
    
    % Epsilon is the probability of choosing one of the less valuable
    % options, reflecting undirected exploration
    p = repmat(epsilon, [1, numel(values)]);
    
    % The best option is choosen with a probability inverse to the amount
    % of undirected exploration
    p(bestoption) = 1 - epsilon;
    
% Softmax function
% ================
elseif strcmpi(fun, 'Softmax')
    
    % Free parameter of the softmax function
    beta = linkspec{2};
    
    % Exponentiate the value of each option weighted by the beta parameter
    ev = exp(beta .* values);
    
    % Compute the sum of the transformed values
    sev = sum(ev);
    
    % Get the probability each choice
    p = ev ./ sev;
end

end