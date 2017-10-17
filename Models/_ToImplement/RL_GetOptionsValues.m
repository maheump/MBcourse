function [ Va, Vb ] = RL_GetOptionsValues( choices, rewards, alpha )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

N = numel(rewards);

rewards = repmat(rewards, [numel(alpha), 1, 2]);
prediction = repmat(0.5, [numel(alpha), 1, 2]);

% For each trial
for t = 2:N
    
    % Compute the prediction error only for the current stim
    currentstim = choices(t);
    otherstim = setdiff([1,2], choices(t));
    delta = rewards(:,t,currentstim) - prediction(:,t-1,currentstim);
    
    % Update the prediction based on this prediction error
    prediction(:,t,currentstim) = prediction(:,t-1,currentstim) + (alpha' .* delta);
    prediction(:,t,otherstim)   = prediction(:,t-1,otherstim);
end

% Need to extract options values
Va = prediction(:,:,1);
Vb = prediction(:,:,2);

end