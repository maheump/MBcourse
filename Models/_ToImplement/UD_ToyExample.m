% Toy example of the "Uncertain Softmax" algorithm
% Maxime Maheu (UNICOG, NeuroSpin, CEA, INSERM & FdV, CRI, Universite Paris Descartes)

clear; close('all');

%% Simulate random data (average p(A) = 0.5)

% Choices and rewards are splitted into runs (each run is a cell in those
% cell arrays)
R = 20; % number of runs
choices = cell(1,R);
rewards = cell(1,R);
for r = 1:R
    N = randi([80,120]); % number of trials in a run
    choices{r} = round(rand(1, N)) + 1; % 1 for A and 2 for B
    rewards{r} = round(rand(1, N)); % 0 for no reward and 1 for a reward
end

%% Define parameters grids

alphagrid   = 0:0.01:2;

epsilongrid = 0:0.01:1;
betagrid    = -10:0.1:10;
beta0grid   = -5:0.1:5;
phigrid     = -10:0.1:10;

%% Run and fit all of the models

tic;
for r = 1:R
disp(r);

% Run the standard Rescorla-Wagner algorithm
[sRL_BIC(r), sRL_LSE(r), sRL_gridsearch{r}, sRL_bestalpha(r), sRL_bestbeta(r), sRL_PofA{r}, sRL_PofB{r}] = ...
    RL_SimulateStandardRL(choices{r}, rewards{r}, alphagrid, betagrid);

% Run the local matching law
[LML_BIC(r), LML_LSE(r), LML_PofA{r}, LML_PofB{r}] = ...
    LML_SimulateLML(choices{r}, rewards{r});

% Run the "epsilon-greedy" choice rule
[EG_BIC(r), EG_LSE(r), EG_gridsearch{r}, EG_bestepsilon(r), EG_PofA{r}, EG_PofB{r}] = ...
    EG_SimulateEG(choices{r}, rewards{r}, epsilongrid);

% Run the classical softmax
[SM_BIC(r), SM_LSE(r), SM_gridsearch{r}, SM_bestbeta(r), SM_PofA{r}, SM_PofB{r}] = ...
    SM_SimulateSoftmax(choices{r}, rewards{r}, betagrid);

% Run the uncertainty-based temperature model
[UBT_BIC(r), UBT_LSE(r), UBT_gridsearch{r}, UBT_bestbeta0(r), UBT_PofA{r}, UBT_PofB{r}] = ...
    UBT_SimulateUBT(choices{r}, rewards{r}, beta0grid);

% Run the uncertain softmax
[US_BIC(r), US_LSE(r), US_gridsearch{r}, US_bestbeta(r), US_bestphi(r), US_PofA{r}, US_PofB{r}] = ...
    US_SimulateUS(choices{r}, rewards{r}, betagrid, phigrid);

end
toc;

%% Mean over runs

sRL_BIC = mean(sRL_BIC);
LML_BIC = nanmean(LML_BIC);
EG_BIC  = mean(EG_BIC);
SM_BIC  = mean(SM_BIC);
UBT_BIC = mean(UBT_BIC);
US_BIC  = mean(US_BIC);

sRL_LSE = mean(sRL_LSE);
LML_LSE = nanmean(LML_LSE);
EG_LSE  = mean(EG_LSE);
SM_LSE  = mean(SM_LSE);
UBT_LSE = mean(UBT_LSE);
US_LSE  = mean(US_LSE);

sRL_bestalpha  = round(mean(sRL_bestalpha));
sRL_bestbeta   = round(mean(sRL_bestbeta));
EG_bestepsilon = round(mean(EG_bestepsilon));
SM_bestbeta    = round(mean(SM_bestbeta));
UBT_bestbeta0  = round(mean(UBT_bestbeta0));
US_bestbeta    = round(mean(US_bestbeta));
US_bestphi     = round(mean(US_bestphi));

sRL_gridsearch = mean(reshape(cell2mat(sRL_gridsearch), [size(cell2mat(sRL_gridsearch),1), size(cell2mat(sRL_gridsearch),2)/R, R]), 3);
EG_gridsearch = mean(reshape(cell2mat(EG_gridsearch), [size(cell2mat(EG_gridsearch),1), size(cell2mat(EG_gridsearch),2)/R, R]), 3);
SM_gridsearch = mean(reshape(cell2mat(SM_gridsearch), [size(cell2mat(SM_gridsearch),1), size(cell2mat(SM_gridsearch),2)/R, R]), 3);
UBT_gridsearch = mean(reshape(cell2mat(UBT_gridsearch), [size(cell2mat(UBT_gridsearch),1), size(cell2mat(UBT_gridsearch),2)/R, R]), 3);
US_gridsearch = mean(reshape(cell2mat(US_gridsearch), [size(cell2mat(US_gridsearch),1), size(cell2mat(US_gridsearch),2)/R, R]), 3);

%% Plot the best parameters

% Default plot parameters
set(0, 'DefaultAxesFontName', 'Times');
set(0, 'DefaultTextFontname', 'Times');
set(0, 'DefaultAxesFontName', 'Times');
set(0, 'DefaultTextFontSize', 15);
set(0, 'DefaultAxesFontSize', 15);
set(0, 'DefaultAxesXColor', 'k');
set(0, 'DefaultAxesYColor', 'k');
set(0, 'DefaultAxesLineWidth', 2);
set(0, 'DefaultAxesLayer', 'Top');
set(0, 'DefaultAxesTickDir', 'out');

% Create a new figure
f1 = figure('Name', 'Quality of fit', 'Units', 'Normalized', 'Position', [0.0229 0.1633 0.9771 0.5925]);
maxLSE = ceil(max([max(sRL_gridsearch(:)); LML_LSE; EG_gridsearch(:); ...
    SM_gridsearch(:); UBT_gridsearch(:); US_gridsearch(:)]))+1;
models = {'Standard RL model', 'Local matching rule', '\epsilon - greedy policy', ...
    'Classical softmax', 'Uncertainty-based temperature model', 'Uncertain softmax'};

% ... for the standard RL model
subplot(2,5,1);
surf(alphagrid, betagrid, sRL_gridsearch, 'EdgeColor', 'None'); hold('on');
plot3(alphagrid(sRL_bestalpha), betagrid(sRL_bestbeta), sRL_LSE, 'o', ...
    'LineWidth', 2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'MarkerSize', 10);
axis([alphagrid(1), alphagrid(end), betagrid(1), betagrid(end), 0, maxLSE]);
view(-70, 35);
xlabel({'', '', '\alpha'}); ylabel('\beta'); zlabel({'Reinforcement-learning models', '', 'SSE'});
title(models{1});

% ... for the local matching law
subplot(2,5,5+1);
bar(1, LML_LSE, 0.25, 'r', 'LineWidth', 2);
ylim([0, maxLSE]); box('off');
set(gca, 'XTick', []); ylabel({'Decision-making models', '', 'SSE'});
title(models{2});

% ... for the epsilon-greedy policy
subplot(2,5,5+2);
plot(epsilongrid, EG_gridsearch, '-', 'LineWidth', 2); hold('on');
plot(epsilongrid(EG_bestepsilon), EG_LSE, 'ko', 'LineWidth', 2, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'MarkerSize', 10);
axis([epsilongrid(1), epsilongrid(end), 0, maxLSE]);
set(gca, 'Box', 'Off');
xlabel('\epsilon'); ylabel('SSE');
title(models{3});

% ... for the classical softmax
subplot(2,5,5+3);
plot(betagrid, SM_gridsearch, '-', 'LineWidth', 2); hold('on');
plot(betagrid(SM_bestbeta), SM_LSE, 'ko', 'LineWidth', 2, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'MarkerSize', 10);
axis([betagrid(1), betagrid(end), 0, maxLSE]);
set(gca, 'Box', 'Off');
xlabel('\beta'); ylabel('SSE');
title(models{4});

% ... for the uncertainty-based temperature model
subplot(2,5,5+4);
plot(beta0grid, UBT_gridsearch, '-', 'LineWidth', 2); hold('on');
plot(beta0grid(UBT_bestbeta0), UBT_LSE, 'ko', 'LineWidth', 2, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'MarkerSize', 10);
axis([beta0grid(1), beta0grid(end), 0, maxLSE]);
set(gca, 'Box', 'Off');
xlabel('\beta_{0}'); ylabel('SSE');
title(models{5});

% ... for the uncertain softmax
subplot(2,5,5+5);
surf(phigrid, betagrid, US_gridsearch, 'EdgeColor', 'None'); hold('on');
plot3(betagrid(US_bestphi), phigrid(US_bestbeta), US_LSE, 'o', ...
    'LineWidth', 2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'MarkerSize', 10);
axis([betagrid(1), betagrid(end), phigrid(1), phigrid(end), 0, maxLSE]);
set(gca, 'Box', 'Off'); view(-110, 38);
xlabel('\phi'); ylabel('\beta'); zlabel('SSE');
title(models{6});

try export_fig('Fig1.png');
catch, end

%% Plot the best model

f2 = figure('Name', 'Best model', 'Units', 'Normalized', 'Position', [0.0865 0.0750 0.8531 0.7700]);
BICs = [sRL_BIC, LML_BIC, EG_BIC, SM_BIC, UBT_BIC, US_BIC];
bar(1:numel(BICs), BICs, 0.25, 'FaceColor', repmat(0.5, 1, 3)); hold('on');
bar(find(BICs == min(BICs)), min(BICs), 0.25, 'r');
set(gca, 'XTick', 1:numel(BICs), 'XTickLabel', models); box('off');
xlabel('Models'); ylabel('BIC'); title('Best model');

try export_fig('Fig2.png');
catch, end
