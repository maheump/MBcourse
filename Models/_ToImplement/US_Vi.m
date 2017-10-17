function Vi = US_Vi( Pi, phi )
%Vi computes the value of an option at each trial given its reward
%probability (Pi) and an uncertainty parameter (phi).
%   "Pi" and "phi" must be 1-row arrays. Their length can be different.

N = numel(Pi);
Pi = repmat(Pi, [numel(phi), 1]);
phi = repmat(phi', [1, N]);

Vi = Pi + (phi .* Pi .* (1-Pi));

end

