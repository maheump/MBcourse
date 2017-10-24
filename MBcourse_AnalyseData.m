% Script to analyze individual data.
% Analysis script for data from a reversal learning task. The data is
% fitted using variants of the RW reinforcement learning models using
% grid-search and likelihood estimation. Different versions of the model
% are compared using standard Bayesian Model Selection techniques.
% 
% Maxime Maheu, 11/2016

%% Initialization
%  ==============

% If the required variables are not in the workspace, load a data file
%datafile = '/Users/Maxime/Documents/MATLAB/MBcourse/SimDataForCourse_StudentsData.mat';
%load(datafile, 'design', 'data', 'volatility');

% Folder in which to save the generated figures
%if ~exist(dirfigsave, 'dir'), dirfigsave = uigetdir; end

%% Display the design and the choices
%  ==================================

% Prepare the window
figure('Color', ones(1,3), 'Units', 'Normalized', 'Position', [0.2328 0.3892 0.5344 0.1825]);
col = mat2cell(winter(6), ones(1,6)); col = col([2,5],:);

% Plot the properties of the design
plot(1:design.nTrials,   design.feedbackprob, '-', 'Color', col{1}, 'LineWidth', 2); hold('on');
plot(1:design.nTrials, 1-design.feedbackprob, '-', 'Color', col{2}, 'LineWidth', 2);
for i = 1:2
    plot(find(((data.choice == i) + (data.choiceOutcome == 1)) == 2), abs(i-3)-1, ...
        'o', 'Color', col{i}, 'MarkerFaceColor', col{i}, 'MarkerSize', 5, 'LineWidth', 1);
    plot(find(((data.choice == i) + (data.choiceOutcome == 0)) == 2), abs(i-3)-1, ...
        'o', 'Color', col{i}, 'MarkerFaceColor', 'w', 'MarkerSize', 5, 'LineWidth', 1);
end

% Customize axes
ylim([-0.1,1.1]);
set(gca, 'Box', 'Off', 'LineWidth', 1, 'FontSize', 15);

% Add some labels
xlabel('Trials');
legend({'p(reward|A)', 'p(reward|B)'}, 'Location', 'EastOutside');

% Save the figure
if exist('volatility', 'var') && ~isempty(volatility)
    ylabel(upper({sprintf('%s volatility', volatility),''}), 'FontWeight', 'Bold');
    save2pdf(fullfile(dirfigsave, sprintf('DesignSimu%s.pdf', upper(volatility))));
else
    ylabel({'YOUR DATA',''}, 'FontWeight', 'Bold');
    save2pdf(fullfile(dirfigsave, 'DesignData.pdf'));
end

%% Model-free analysis: performance
%  ================================

% Get averaged subject's performance
M = mean(data.choiceOutcome);
S = std(data.choiceOutcome) / sqrt(numel(data.choiceOutcome));

% Prepare the window
figure('Color', ones(1,3), 'Units', 'Normalized', 'Position', [0.4479 0.2508 0.1047 0.4600]);

% Plot the results
bar(M, 'LineWidth', 1, 'FaceColor', [0.2078 0.1647 0.5255]); hold('on');
plot(ones(1,2), [M-S, M+S], 'k-', 'LineWidth', 1);
axis([0.5,1.5,0,1]);
set(gca, 'XTick', [], 'Box', 'Off', 'LineWidth', 1, 'FontSize', 15);
ylabel('Averaged performance (+/- SEM)');

% Save the figure
if exist('volatility', 'var') && ~isempty(volatility)
    title(upper({sprintf('%s volatility', volatility),''}));
    save2pdf(fullfile(dirfigsave, sprintf('PerfSimu%s.pdf', upper(volatility))));
else
    title({'YOUR DATA',''});
    save2pdf(fullfile(dirfigsave, 'PerfData.pdf'));
end

% % Temptative of reversal evoked performance (not enough reversals though)
% rev = find(diff(design.feedbackprob) ~= 0);
% Nrev = numel(rev);
% revperf = NaN(Nrev,10+20+1);
% for i = 1:Nrev
%     revperf(i,:) = data.outcome(rev(i)-10:rev(i)+20);
% end

%% Model-based analysis: grid-search
%  =================================

% Define grids
reso = 100;
alpha   = linspace(0.1, 0.9, reso);
beta    = linspace(0.1, 8,   reso);
epsilon = linspace(0.1, 0.9, reso);

% Prepare outputs
MSE = cell(1,3);
LLH = cell(1,3);
LH  = cell(1,3);

% For each parameters set
for a = 1:reso
    for b = 1:reso
        
        % Simulate RL model
        SimuSoftmax = MBcourse_RLobs_Simulation(alpha(a), {'Softmax', beta(b)},   design.feedback, 1, data.choice);
        SimuGreedy  = MBcourse_RLobs_Simulation(alpha(a), {'Greedy', epsilon(b)}, design.feedback, 1, data.choice);
        SimuLMR     = MBcourse_RLobs_Simulation(alpha(a), {'LMR'},                design.feedback, 1, data.choice);

        % Measure quality of fit
        [MSE{1}(b,a), LLH{1}(b,a)] = MBcourse_RLobs_Fit(SimuSoftmax.choiceProba, data.choice);
        [MSE{2}(b,a), LLH{2}(b,a)] = MBcourse_RLobs_Fit(SimuGreedy.choiceProba, data.choice);
    end
	[MSE{3}(a), LLH{3}(a)] = MBcourse_RLobs_Fit(SimuLMR.choiceProba, data.choice);
end

% Compute the delta-likelihood (in the log-space)
for i = 1:3, LH{i} = exp(LLH{i}); end

%% Grid-search with LLH
%  ====================

% Prepare the figure
figure('Color', ones(1,3), 'Units', 'Normalized', 'Position', [0.2547 0.3842 0.4911 0.1933]);

% Display the log-likelihood
subplot(1,3,1);
imagesc(alpha, beta, LLH{1}); hold('on');
contour(alpha, beta, LLH{1}, 'k-');
cbr = colorbar; cbr.LineWidth = 1; cbr.Label.FontSize = 15;
axis('square'); axis('square'); axis('xy');
set(gca, 'FontSize', 15);
xlabel('Grid for the \alpha parameter');
ylabel('Grid for the \beta parameter');
title('Log-likelihood');

% Find the best parameters
[~, bestidx] = max(LLH{1}(:));
[best_b, best_a] = ind2sub(size(LH{1}), bestidx);
best_alpha = alpha(best_a);
best_beta = beta(best_b);
plot(best_alpha, best_beta, 'k.', 'MarkerSize', 12);
text(best_alpha, best_beta, ['$\hat{\alpha} = ', num2str(best_alpha, 2), ...
    ', \hat{\beta} = ', num2str(best_beta, 2), '$'], 'Interpreter', 'LaTeX', ...
    'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Bottom');

% Display the likelihood
subplot(1,3,2);
imagesc(alpha, beta, LH{1}); hold('on');
contour(alpha, beta, LH{1}, 'k-');
cbr = colorbar; cbr.LineWidth = 1; cbr.Label.FontSize = 15;
axis('square'); axis('square'); axis('xy');
set(gca, 'FontSize', 15);
xlabel('Grid for the \alpha parameter');
ylabel('Grid for the \beta parameter');
title('Likelihood');

% Find the best parameters
[~, bestidx] = max(LH{1}(:));
[best_b, best_a] = ind2sub(size(LH{1}), bestidx);
best_alpha = alpha(best_a);
best_beta = beta(best_b);
plot(best_alpha, best_beta, 'k.', 'MarkerSize', 12);
text(best_alpha, best_beta, ['$\hat{\alpha} = ', num2str(best_alpha, 2), ...
    ', \hat{\beta} = ', num2str(best_beta, 2), '$'], 'Interpreter', 'LaTeX', ...
    'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Bottom');

% Display the mean-squared error
subplot(1,3,3);
imagesc(alpha, beta, MSE{1}); hold('on');
contour(alpha, beta, MSE{1}, 'k-');
cbr = colorbar; cbr.LineWidth = 1; cbr.Label.FontSize = 15;
axis('square'); axis('square'); axis('xy');
set(gca, 'FontSize', 15);
xlabel('Grid for the \alpha parameter');
ylabel('Grid for the \beta parameter');
title('Mean squared error');

[~, bestidx] = min(MSE{1}(:));
[best_b, best_a] = ind2sub(size(LH{1}), bestidx);
best_alpha = alpha(best_a);
best_beta = beta(best_b);
plot(best_alpha, best_beta, 'k.', 'MarkerSize', 12);
text(best_alpha, best_beta, ['$\hat{\alpha} = ', num2str(best_alpha, 2), ...
    ', \hat{\beta} = ', num2str(best_beta, 2), '$'], 'Interpreter', 'LaTeX', ...
    'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Bottom');

% Customize the plot
axes('Units', 'normalized', 'Position', [0 0 1 1])
if exist('volatility', 'var') && ~isempty(volatility)
    txt = upper({sprintf('%s\n volatility', volatility),''});
else
    txt = 'YOUR DATA';
end
text(-0.9, 0, txt, 'FontSize', 20, 'FontWeight', 'Bold', 'Rotation', 90, ...
    'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
axis('off'); axis([-1,1,-1,1]);

% Save the figure
if exist('volatility', 'var') && ~isempty(volatility)
    save2pdf(fullfile(dirfigsave, sprintf('GsSimu%s.pdf', upper(volatility))));
else
    save2pdf(fullfile(dirfigsave, 'GsData.pdf'));
end

%% Compare models
%  ==============

% Prepare the figure
figure('Color', ones(1,3), 'Units', 'Normalized', 'Position', [0.2266 0.3783 0.5469 0.2042]);
models = {'RL + Softmax', 'RL + \epsilon-greedy', 'RL + LMR'};

% Display the likelihood of the first model
subplot(1,4,1);
imagesc(alpha, beta, LH{1});
axis('square'); axis('xy');
set(gca, 'FontSize', 15, 'LineWidth', 1);
xlabel('\alpha'); ylabel('\beta');
title({models{1},''});

% Display the likelihood of the second model
subplot(1,4,2);
fill([alpha(1), alpha, alpha(end)], [0, LH{3}, 0], 'k', 'EdgeColor', lines(1), 'FaceColor', lines(1), 'LineWidth', 2, 'FaceAlpha', 0.3);
axis('square');
set(gca, 'FontSize', 15, 'LineWidth', 1);
xlabel('\alpha'); ylabel('Likelihood');
title({models{3},''});

% Compare models based on maximum likelihood estimates
subplot(1,4,3);
L = NaN(1,2);
L(1) = max(log(LH{1}(:)));
bar(1, L(1), 'FaceColor', repmat(0.6,1,3), 'LineWidth', 1); hold('on');
L(2) = max(log(LH{3}(:)));
bar(2, L(2), 'FaceColor', repmat(0.6,1,3), 'LineWidth', 1); hold('on');
[L,i] = nanmax(L);
bar(i, L, 'FaceColor', 'r', 'LineWidth', 1); hold('on');
axis('square'); xlim([0,3]);
set(gca, 'XTick', [1,2], 'XTickLabel', models([1,3]), 'XTickLabelRotation', 20);
set(gca, 'FontSize', 15, 'LineWidth', 1);
ylabel('Log-likelihood'); title({'Maximum likelihood', '(accuracy)'});

% Compare models based on models' evidence
subplot(1,4,4);
L = NaN(1,2);
L(1) = mean(log(LH{1}(:)));
bar(1, L(1), 'FaceColor', repmat(0.6,1,3), 'LineWidth', 1); hold('on');
L(2) = mean(log(LH{3}(:)));
bar(2, L(2), 'FaceColor', repmat(0.6,1,3), 'LineWidth', 1); hold('on');
[L,i] = nanmax(L);
bar(i, L, 'FaceColor', 'r', 'LineWidth', 1); hold('on');
axis('square'); xlim([0,3]);
set(gca, 'XTick', [1,2], 'XTickLabel', models([1,3]), 'XTickLabelRotation', 20);
set(gca, 'FontSize', 15, 'LineWidth', 1);
ylabel('Log-likelihood'); title({'Model evidence','(accuracy + complexity)'});

% Customize the plot
axes('Units', 'normalized', 'Position', [0 0 1 1])
if exist('volatility', 'var') && ~isempty(volatility)
    txt = upper({sprintf('%s\n volatility', volatility),''});
else
    txt = 'YOUR DATA';
end
text(-0.9, 0, txt, 'FontSize', 20, 'FontWeight', 'Bold', 'Rotation', 90, ...
    'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
axis('off'); axis([-1,1,-1,1]);

% Save the figure
if exist('volatility', 'var') && ~isempty(volatility)
    save2pdf(fullfile(dirfigsave, sprintf('BmsSimu%s.pdf', upper(volatility))));
else
    save2pdf(fullfile(dirfigsave, 'BmsData.pdf'));
end

%% Run simulations with the best parameters
%  ========================================

% Simulate a RL observer with these parameters
simu = MBcourse_RLobs_Simulation(best_alpha, {'Softmax', best_beta}, design.feedback, 1, data.choice);

% Prepare the window
figure('Color', ones(1,3), 'Units', 'Normalized', 'Position', [0.2328 0.3892 0.5344 0.1825]);
col = mat2cell(winter(6), ones(1,6)); col = col([2,5],:);

% Plot the properties of the design
plot(1:design.nTrials,   design.feedbackprob, '--', 'Color', col{1}, 'LineWidth', 1); hold('on');
plot(1:design.nTrials, 1-design.feedbackprob, '--', 'Color', col{2}, 'LineWidth', 1);

% Plot stimuli values
plot(1:design.nTrials, simu.values(:,1), '-', 'Color', col{1}, 'LineWidth', 2);
plot(1:design.nTrials, simu.values(:,2), '-', 'Color', col{2}, 'LineWidth', 2);

% Plot the choice probability
for i = 1:2, plot(1:design.nTrials, simu.choiceProba(:,i), '.', 'Color', col{i}, 'MarkerSize', 10, 'LineWidth', 1); end

% Plot choices
for i = 1:2
    plot(find(((simu.choice == i) + (simu.choiceOutcome == 1)) == 2), abs(i-3)-1, ...
        'o', 'Color', col{i}, 'MarkerFaceColor', col{i}, 'MarkerSize', 5, 'LineWidth', 1);
    plot(find(((simu.choice == i) + (simu.choiceOutcome == 0)) == 2), abs(i-3)-1, ...
        'o', 'Color', col{i}, 'MarkerFaceColor', 'w', 'MarkerSize', 5, 'LineWidth', 1);
end

% Customize the axes
ylim([-0.1, 1.1]);
set(gca, 'FontSize', 15, 'LineWidth', 1);

% Add some labels
legend({'p(reward|A)', 'p(reward|B)', 'Value of A', 'Value of B', 'p(choose A)', 'p(choose B)'}, 'Location', 'EastOutside');
xlabel('Trials'); ylabel('Probability');

% Save the figure
if exist('volatility', 'var') && ~isempty(volatility)
    ylabel(upper({sprintf('%s volatility', volatility),''}), 'FontWeight', 'Bold');
    save2pdf(fullfile(dirfigsave, sprintf('BestSimu%s.pdf', upper(volatility))));
else
    ylabel({'YOUR DATA',''}, 'FontWeight', 'Bold');
    save2pdf(fullfile(dirfigsave, 'BestSimuData.pdf'));
end

%% Compute prediction error
%  =========================

% Prepare the window
figure('Color', ones(1,3), 'Units', 'Normalized', 'Position', [0.2328 0.3892 0.5344 0.1825]);

% Map depicting positive/negative PE
x = linspace(-1, 1, reso);
imagesc(1:design.nTrials, x, repmat(x', 1, design.nTrials), 'AlphaData', 1/2); hold('on');
plot([1, design.nTrials], zeros(1,2), 'k--', 'LineWidth', 1);
colormap(cbrewer2('RdBu'));

% Display prediction error levels from the observer
plot(1:design.nTrials, simu.predictionErrors, 'k-', 'LineWidth', 3);

% Customize the axes
axis('xy');
set(gca, 'FontSize', 15, 'LineWidth', 1);

% Add some labels
xlabel('Trials'); ylabel('Prediction error');

if exist('volatility', 'var') && ~isempty(volatility)
    ylabel(upper({sprintf('%s volatility', volatility),''}), 'FontWeight', 'Bold');
    save2pdf(fullfile(dirfigsave, sprintf('Pe%s.pdf', upper(volatility))));
else
    ylabel({'YOUR DATA',''}, 'FontWeight', 'Bold');
    save2pdf(fullfile(dirfigsave, 'PeData.pdf'));
end
