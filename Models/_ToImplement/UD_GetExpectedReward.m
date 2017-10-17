function [ pa, pb, N ] = UD_GetExpectedReward( choices, rewards )
%UD_GetExpectedReward computes expected reward rates for binary choices
%based on their related reward probabilities.
%	In the US, where the choices are at steady-state, the expected reward
%	is taken as the reward probability. We just have to compute it
%	iteratively using rewards.

N = numel(choices);

pa = NaN(1,N);
pb = NaN(1,N);

for t = 1:N
    pa(t) = mean(rewards(choices(1:t) == 1));
    pb(t) = mean(rewards(choices(1:t) == 2));
end

pa(isnan(pa)) = 0.5;
pb(isnan(pb)) = 0.5;

end