% Routine to run a reversal learning experiment, save, preprocess, analyse
% the data and update the slides of the course with the newly obtained
% figures.
% 
% Maxime Maheu, October 2016

%% SETUP PSYCHTOOLBOX
%  ==================

% Clear the place
clear; close('all'); clc;

% Define the MATLAB directory
matlabrepo = '/Users/Maxime/Documents/MATLAB';
cd(fullfile(matlabrepo, 'Psychtoolbox'));
SetupPsychtoolbox;
fprintf('\n=====> PSYCHTOOLBOX IS IN THE PATH\n');

%% RUN THE EXPERIMENT
%  ==================

% Make sure the functions are in the MATLAB path
expedir = fullfile(matlabrepo, 'RLtutorial_codeNdata', 'taskReversalLearn');
cd(expedir); addpath(genpath(pwd));

% Run the experiment
% - odd subject numbers: high volatility (volatile) condition
% - even subject numbers: low volatility (stable) condition.  
s = 50;
revlRun(s);
fprintf('\n=====> DATA ACQUIRED\n');

%% MOVE THE DATA FILE TO THE COURSE 
%  ================================

% Locate the data file
datdir = fullfile(expedir, 'data');
datafile = dir(fullfile(datdir, sprintf('tutorialRevLearn_*_s%03.0f_data*.mat', s)));
datafile = fullfile(datdir, datafile.name);

% Define the 
coursedir = fullfile(matlabrepo, 'MBcourse');
datafolder = 'Data_2017';
mkdir(fullfile(coursedir, datafolder));

% Copy the file to the 
datelab = datestr(now, 'ddmmyy');
mvdatafile = fullfile(coursedir, datafolder, sprintf('DataForCourse_%s.mat', datelab));
[success, msg, msgid] = copyfile(datafile, mvdatafile);
fprintf('\n=====> DATA FILE MOVED TO ANALYSIS REPOSITORY\n');

%% PREPROCESS THE DATA
%  ===================

% Load the data
cd(coursedir); addpath(genpath(pwd));
load(mvdatafile);

% Create the design structure
design.nTrialsPerBlocks = diff([0;find(diff(data.prep.feedbackprob)~=0);data.prep.nt])';
design.pReward          = data.prep.prob;
design.nTrials          = data.prep.nt;
design.nBlocks          = sum(diff([data.prep.feedbackprob;1]) ~= 0);
design.probTrial        = data.prep.feedbackprob([find(diff(data.prep.feedbackprob) ~= 0);data.prep.nt])';
design.feedbackprob     = data.prep.feedbackprob;
design.feedback         = data.prep.feedback;

% Create the data structure
data.choice        = data.choice;
data.choiceOutcome = data.outcome;

% This variable is empty to make sure that other scripts can work
volatility = [];

% Save these 3 variables
ppdatafile = fullfile(coursedir, datafolder, 'SimDataForCourse_StudentsData.mat');
save(ppdatafile, 'data', 'design', 'volatility');
fprintf('\n=====> DATA SUCCESSFULY PREPROCESSED\n');

%% RUN THE ANALYSIS SCRIPT
%  =======================

% Specify where to save 
coursedir2 = '/Users/Maxime/Documents/My courses/A hands-on introduction to the model-based approach in neuroscience/';
dirfigsave = fullfile(coursedir2, 'data');

% Perform model-free and model-based analyses
MBcourse_AnalyseData;
close('all');

%% UPDATE THE SLIDES WITH THE FIGURES ABOUT THE NEW DATAFILE
%  =========================================================

% Add pdflatex functions to the path
setenv('PATH', [getenv('PATH') ':/Library/TeX/texbin']);

% Replace spaces " " by "\ " to ensure LaTeX is working
coursedir2 = strrep(coursedir2, ' ', '\ ');

% Update the slides with the new figures using LaTeX composition
slidesfile = 'Computational_neuroscience_course.tex';
cmd = sprintf('cd %s && pdflatex %s', coursedir2, slidesfile);
status = system(cmd, '-echo');
if ~status, fprintf('\n=====> THE SLIDES WERE SUCCESSFULY UPDATED\n'); end

% Open the updated slides
disp('=====> THE SLIDES ARE ABOUT TO BE DISPLAYED');
cmd = sprintf('open %s', fullfile(coursedir2, [slidesfile(1:end-3) 'pdf']));
[status, result] = system(cmd, '-echo');
