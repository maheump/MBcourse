%% Simulated data
%  ==============

% Low-volatility experiment
design = MBcourse_GenerateTaskDesign([70 70], 0.7);
data = MBcourse_RLobs_Simulation(0.55, {'Softmax', 4}, design.feedback, 1);
volatility = 'low';
save('SimDataForCourse_LowVolatility.mat', 'data', 'design', 'volatility');

% High-volatility experiment
design = MBcourse_GenerateTaskDesign([40 20 30 50], 0.7);
data = MBcourse_RLobs_Simulation(0.55, {'Softmax', 4}, design.feedback, 1);
volatility = 'high';
save('SimDataForCourse_HighVolatility.mat', 'data', 'design', 'volatility');

%% Students' data
%  ==============

% Load the data
d = load('DataForCourse_191116.mat');

% Create the design structure
design.nTrialsPerBlocks = diff([0;find(diff(d.data.prep.feedbackprob)~=0);d.data.prep.nt])';
design.pReward          = d.data.prep.prob;
design.nTrials          = d.data.prep.nt;
design.nBlocks          = sum(diff([d.data.prep.feedbackprob;1]) ~= 0);
design.probTrial        = d.data.prep.feedbackprob([find(diff(d.data.prep.feedbackprob) ~= 0);d.data.prep.nt])';
design.feedbackprob     = d.data.prep.feedbackprob;
design.feedback         = d.data.prep.feedback;

% Create the data structure
data.choice        = d.data.choice;
data.choiceOutcome = d.data.outcome;

% This variable is empty to make sure that other scripts can work
volatility = [];

% Save these 3 variables
save('SimDataForCourse_StudentsData.mat', 'data', 'design', 'volatility');