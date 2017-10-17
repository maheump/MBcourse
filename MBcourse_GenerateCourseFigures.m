% Script generating the figures for the course.
% 
% Maxime Maheu, 11/2016

%% Investigating the alpha parameter (1)
%  =====================================

% Get task design
N = 50;
fb_acq = [ones(N,1); zeros(N,1)];
fb_acq = [fb_acq, fb_acq];

% Define parameters with which simulate the agents
alpha = 0.1:0.05:0.9;

% Prepare outputs
nA = numel(alpha);
simu = cell(nA, 1);

% For each value of alpha parameter
for a = 1:nA
    
    % Simulate the agents
    simu{a} = MBcourse_RLobs_Simulation(alpha(a), {'LMR'}, fb_acq, 1000);
    
    % Average over agents
    simu{a}.choiceOutcome = mean(simu{a}.choiceOutcome, 2);
    simu{a}.values        = mean(simu{a}.values, 3);
    simu{a}.choiceProba   = mean(simu{a}.choiceProba, 3);
end

% Prepare figure
figure('Position', [0.2703 0.3400 0.4599 0.2808]);
cols = flipud(autumn(nA));

% Acquisition
subplot(1,2,1);
plot([1,N], ones(1,2)./2, 'k--', 'LineWidth', 1); hold('on');
for a = 1:nA
    plot(1:N, simu{a}.values(1:N,1), '.-', 'Color', cols(a,:), ...
        'MarkerSize', 12, 'LineWidth', 2);
end
axis([1,N,0,1]); axis('square'); grid('on');
colormap(flipud(autumn)); cbr = colorbar; caxis([min(alpha), max(alpha)]);
cbr.LineWidth = 1; cbr.Label.String = 'Value of \alpha parameter'; cbr.Label.FontSize = 15;
set(gca, 'FontSize', 15, 'LineWidth', 1, 'Layer', 'Bottom');
xlabel('Trials'); ylabel('Value of option A'); title({'Acquisition phase','(value of A = 1)'});

% Extinction
subplot(1,2,2);
plot([(N+1),(2*N)], ones(1,2)./2, 'k--', 'LineWidth', 1); hold('on');
for a = 1:nA
    plot((N+1):(2*N), simu{a}.values((N+1):(2*N),1), '.-', 'Color', cols(a,:), ...
        'MarkerSize', 12, 'LineWidth', 2); hold('on');
end
colormap(flipud(autumn)); cbr = colorbar; caxis([min(alpha), max(alpha)]);
cbr.LineWidth = 1; cbr.Label.String = 'Value of \alpha parameter'; cbr.Label.FontSize = 15;
axis([(N+1),(2*N),0,1]); axis('square'); grid('on');
set(gca, 'FontSize', 15, 'LineWidth', 1, 'Layer', 'Bottom');
xlabel('Trials'); ylabel('Value of option A'); title({'Extinction phase','(value of A = 0)'});

%% Investigating the alpha parameter (2)
%  =====================================

% Get task-design
design = MBcourse_GenerateTaskDesign([30 10 20 40], 0.8);

% Define parameters with which simulate the agents
alpha = [0.1, 0.2, 0.3, 0.5, 0.7];

% Prepare outputs
nA = numel(alpha);
simu = cell(nA, 1);

% For each value of alpha parameter
for a = 1:nA
    
    % Simulate the agents
    simu{a} = MBcourse_RLobs_Simulation(alpha(a), {'LMR'}, design.feedback, 1000);
    
    % Average over agents
    simu{a}.choiceOutcome = mean(simu{a}.choiceOutcome, 2);
    simu{a}.values        = mean(simu{a}.values, 3);
    simu{a}.choiceProba   = mean(simu{a}.choiceProba, 3);
end

% Prepare figure
figure('Position', [0.2635 0.3558 0.4729 0.2492]);
cols = flipud(autumn(nA));

% Plot models' choices
plot([1,design.nTrials], ones(1,2)./2, 'k--', 'LineWidth', 1); hold('on');
plot(1:design.nTrials, design.feedbackprob, 'k-', 'LineWidth', 1);
for a = 1:nA
    plot(1:design.nTrials, simu{a}.values(:,1), '.-', 'Color', cols(a,:), ...
        'MarkerSize', 12, 'LineWidth', 2);
end

% Customize the axes
colormap(flipud(autumn)); cbr = colorbar; caxis([min(alpha), max(alpha)]);
cbr.LineWidth = 1; cbr.Label.FontSize = 15;
axis([1,design.nTrials,0,1]); grid('on');
set(gca, 'FontSize', 15, 'LineWidth', 1, 'Layer', 'Bottom');

% Add some labels
cbr.Label.String = 'Value of \alpha parameter';
xlabel('Trials'); ylabel('Value of option A'); title({'Volatile environment',''});

%% Investigating the beta parameter
%  ================================
 
% Choose a few values for the beta parameter
beta = [3 5 10 20];
nB = numel(beta);

% Grid for values
vA = 0:0.01:1;
vB = 1 - vA;

% Useful variables
ct = 0;
x = .04;

% Prepare output variables
lgd1 = cell(1, nB);
lgd2 = NaN(1, nB);
pA = NaN(nB, length(vA));
choices = NaN(nB, length(vA));

% Loop over possible values of beta
for beta = beta
    ct = ct+1;
    
    % Compute p(choosing A)
    pA(ct,:) = exp(beta*vA) ./ (exp(beta*vA) + exp(beta*vB));
    lgd1{ct} = ['\beta = ', sprintf('%d',beta)];
    
    % Generate choices based on this data
    tmp = double(rand(1,length(vA)) < pA(ct,:))*(1+x*ct);
    tmp(tmp == 0) = tmp(tmp == 0) - x * ct;
    choices(ct,:) = tmp;
end

% Prepare the figure
figure('Position', [0.2745 0.3617 0.4516 0.2383]);
cols = lines(nB);

% Plot helping lines
plot([-1, 1], ones(1,2)./2, 'k--', 'LineWidth', 1); hold('on');
plot(zeros(1,2), [-x*ct, 1+x*ct], 'k--', 'LineWidth', 1);

% Plot the choices
for b = 1:nB
    lgd2(b) = plot(vA-vB, pA(b,:), 'Color', cols(b,:), 'LineWidth', 2);
    plot(vA-vB, choices(b,:), 'ko', 'MarkerFaceColor', cols(b,:), 'LineWidth', 1);
end

% Customize the axes
ylim([-x*ct, 1+x*ct]); grid('on');
set(gca, 'FontSize', 15, 'LineWidth', 1, 'Layer', 'Bottom');

% Add some labels
legend(lgd2, lgd1, 'Location', 'East');
xlabel('V(A) - V(B)');
ylabel('p(choice = A)');

%% Display the influence of alpha parameter on RL-based learning
%  ==============================================================

% Get task design
design = MBcourse_GenerateTaskDesign([30 10 20 40], 0.7);

% Number of subjects to simulate
nSimu = 2000;

% Define parameters with which simulate the agents
alpha = [.1, .4, 1];
beta  = 9;

% Prepare outputs
nA = numel(alpha);
nB = numel(beta);
simu = cell(nA, nB);

% For each parameters set
for a = 1:nA
    for b = 1:nB
        
        % Simulate the agents
        simu{a,b} = MBcourse_RLobs_Simulation(alpha(a), {'Softmax', beta(b)}, design.feedback, nSimu);
        
        % Average over agents
        simu{a,b}.choiceOutcome = mean(simu{a,b}.choiceOutcome, 2);
        simu{a,b}.values        = mean(simu{a,b}.values, 3);
        simu{a,b}.choiceProba   = mean(simu{a,b}.choiceProba, 3);
    end
end

% Plot the results
figure;
for a = 1:nA
    for b = 1:nB
        subplot(nA,nB,b+(nB*(a-1)));
        plot(1:design.nTrials, simu{a,b}.choiceOutcome, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 3); hold('on');
        plot(1:design.nTrials, design.feedbackprob, 'k-');
        plot(1:design.nTrials, simu{a,b}.values, '--', 'LineWidth', 1);
        plot(1:design.nTrials, simu{a,b}.choiceProba(:,1), 'LineWidth', 1);
        
        % Customize the axes
        set(gca, 'FontSize', 15, 'LineWidth', 1);
        axis([1, design.nTrials, 0, 1]);
        
        % Add some labels
        if a == nA && b == nB, legend({'Choice', 'p(reward|A)', 'V(A)', 'V(B)', 'p(choose A)'}); end
        if a == nA, xlabel('Trials'); end
        ylabel('Probability');
        title(['\alpha', sprintf(' = %1.2f', alpha(a)), ', \beta', ...
            sprintf(' = %1.2f', beta(b))], 'FontWeight', 'Normal');
    end
end

%% Plot simulations of one particular agent
%  ========================================

% Get task design
design = MBcourse_GenerateTaskDesign([40 20 30 50], 0.7);

% Simulate the agent
simu = MBcourse_RLobs_Simulation(0.25, {'Softmax', 4}, design.feedback, 1);

% Average over agents
simu.choiceOutcome = mean(simu.choiceOutcome, 2);
simu.values        = mean(simu.values, 3);
simu.choiceProba   = mean(simu.choiceProba, 3);

% Plot the simulations
figure('Position', [0.2714 0.3808 0.4573 0.2000]);
plot(1:design.nTrials, simu.choiceOutcome, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 3); hold('on');
plot(1:design.nTrials, design.feedbackprob, 'k-');
plot(1:design.nTrials, simu.choiceProba(:,1), 'b-', 'LineWidth', 2);
plot(1:design.nTrials, simu.values(:,1), 'c-', 'LineWidth', 0.5);
plot(1:design.nTrials, simu.values(:,2), 'r-', 'LineWidth', 0.5);

% Customize the axes
set(gca, 'FontSize', 15, 'LineWidth', 1);

% Add some labels
legend({'Choice', 'p(reward|A)', 'p(choose A)', 'V(A)', 'V(B)'}, 'Location', 'EastOutside');
xlabel('Trials'); ylabel('Probability');

%% Plot simulations of a group of agents
%  =====================================

% Get task design
design = MBcourse_GenerateTaskDesign([40 20 30 50], 0.7);

% Simulate the agent
simu = MBcourse_RLobs_Simulation(0.25, {'Softmax', 4}, design.feedback, 100);

% Average over agents
simu.choice        = mean(simu.choice, 2);
simu.choiceOutcome = mean(simu.choiceOutcome, 2);
simu.values        = mean(simu.values, 3);
simu.choiceProba   = mean(simu.choiceProba, 3);

% Plot the simulations
figure('Position', [0.2714 0.3808 0.4573 0.2000]);
plot(1:design.nTrials, design.feedbackprob, 'k-'); hold('on');
plot(1:design.nTrials, abs(simu.choice-2), 'k-', 'LineWidth', 2);
plot(1:design.nTrials, simu.choiceProba(:,1), 'b-', 'LineWidth', 2);
plot(1:design.nTrials, simu.values(:,1), 'c-', 'LineWidth', 0.5);
plot(1:design.nTrials, simu.values(:,2), 'r-', 'LineWidth', 0.5);

% Customize the axes
set(gca, 'FontSize', 15, 'LineWidth', 1);

% Add some labels
legend({'p(reward|A)', 'Choice', 'p(choose A)', 'V(A)', 'V(B)'}, 'Location', 'EastOutside');
xlabel('Trials'); ylabel('Probability');
