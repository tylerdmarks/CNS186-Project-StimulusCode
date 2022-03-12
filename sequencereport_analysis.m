

%% import data
clear
fns = uigetfile('.mat', 'MultiSelect', 'on');
num_subjects = length(fns);
for ff = 1:num_subjects
    data(ff) = importdata(fns{ff});
end

% outlier statistics
rt_thresh = 10;         % if first response time is greater than 10 seconds, exclude the trial

%% Response time
f0 = figure;
f1 = figure;
offset = 0.2;
pooled_ft_inf = [];
pooled_ft_sup = [];
pooled_tt_inf = [];
pooled_tt_sup = [];
for ii = 1:num_subjects
    expt = data(ii).expt;
    response = data(ii).response;

    ft_inf = response.firsttime(expt.field == 0);
    ft_sup = response.firsttime(expt.field == 1);
    tt_inf = response.totaltime(expt.field == 0);
    tt_sup = response.totaltime(expt.field == 1);
    
    % exclude outliers
    ft_inf = ft_inf(ft_inf <= rt_thresh);
    ft_sup = ft_sup(ft_sup <= rt_thresh);
    tt_inf = tt_inf(tt_inf <= rt_thresh + 5);
    tt_sup = tt_sup(tt_sup <= rt_thresh + 5); 
    
    figure(f0)
    violinplot(ii-offset, ft_inf', {sprintf('%d', ii)}, 'Width', 0.2, 'ViolinColor', [0.5 0.8 0.1], 'DataAlpha', 0.4, 'ViolinAlpha', 0.1);
    violinplot(ii+offset, ft_sup', {sprintf('%d', ii)}, 'Width', 0.2, 'ViolinColor', [0.3 0.2 0.9], 'DataAlpha', 0.4, 'ViolinAlpha', 0.1);
    p_firsttime(ii) = ranksum(ft_inf, ft_sup);
    text(ii, max([max(ft_inf) max(ft_sup)]), sprintf('p = %0.3f', p_firsttime(ii)));
    title('First response time')
      
    figure(f1)
    violinplot(ii-offset, tt_inf', {sprintf('%d', ii)}, 'Width', 0.2, 'ViolinColor', [0.5 0.8 0.1], 'DataAlpha', 0.4, 'ViolinAlpha', 0.1);
    violinplot(ii+offset, tt_sup', {sprintf('%d', ii)}, 'Width', 0.2, 'ViolinColor', [0.3 0.2 0.9], 'DataAlpha', 0.4, 'ViolinAlpha', 0.1);   
    p_totaltime(ii) = ranksum(tt_inf, tt_sup);
    text(ii, max([max(tt_inf) max(tt_sup)]), sprintf('p = %0.3f', p_totaltime(ii)));
    title('Total response time')
    
    % pool average data across individuals
    pooled_ft_inf = [pooled_ft_inf mean(ft_inf)];
    pooled_ft_sup = [pooled_ft_sup mean(ft_sup)];
    pooled_tt_inf = [pooled_tt_inf mean(tt_inf)];
    pooled_tt_sup = [pooled_tt_sup mean(tt_sup)];  
end
figure(f0)
xlim([0.5 num_subjects+0.5])
xlabel('Participant')
ylabel('First response time (s)')
% legend(plotdata, {'no primer', 'primer'})
title('Sequence task, by participant')
xticks(1:num_subjects)
for ss = 1:num_subjects
    labels{ss} = num2str(ss);
end
xticklabels(labels);

figure(f1)
xlim([0.5 num_subjects+0.5])
xlabel('Participant')
ylabel('Total response time (s)')
% legend(plotdata, {'no primer', 'primer'})
title('Sequence task, by participant')
xticks(1:num_subjects)
for ss = 1:num_subjects
    labels{ss} = num2str(ss);
end
xticklabels(labels);


%%%% Plotting  pooled data

% firsttime
ftinf_avg = mean(pooled_ft_inf);
ftinf_se = std(pooled_ft_inf)/sqrt(length(pooled_ft_inf));
ftsup_avg = mean(pooled_ft_sup);
ftsup_se = std(pooled_ft_sup)/sqrt(length(pooled_ft_sup));

bar_vector = [ftinf_avg ftsup_avg];
error_vector = [ftinf_se ftsup_se];
[~, pooled_p] = ttest(pooled_ft_inf, pooled_ft_sup);

figure
hold on
% bar(1:2, bar_vector, 1);
for dd = 1:length(pooled_ft_inf)
    plot(1:2, [pooled_ft_inf(dd) pooled_ft_sup(dd)], 'k');
    scatter(1, pooled_ft_inf(dd), 'filled', 'MarkerFaceColor', 'k');
    scatter(2, pooled_ft_sup(dd), 'filled', 'MarkerFaceColor', 'k');
end
scatter(1, ftinf_avg, 'filled', 'r');
scatter(2, ftsup_avg, 'filled', 'r');
errorbar(1:2, bar_vector, error_vector, '.');
text(1.5, max([ftinf_avg ftsup_avg]) + 0.025, sprintf('p = %0.3f', pooled_p));
xlabel('Condition')
xlim([0 3])
ylim([0 3])
ylabel('Mean first response time (s)')
xticks([1 2])
xticklabels({'Inf', 'Sup'})
title('Sequence task, pooled data')
axis square

% totaltime
ttinf_avg = mean(pooled_tt_inf);
ttinf_se = std(pooled_tt_inf)/sqrt(length(pooled_tt_inf));
ttsup_avg = mean(pooled_tt_sup);
ttsup_se = std(pooled_tt_sup)/sqrt(length(pooled_tt_sup));

bar_vector = [ttinf_avg ttsup_avg];
error_vector = [ttinf_se ttsup_se];
[~, pooled_p] = ttest(pooled_tt_inf, pooled_tt_sup);

figure
hold on
% bar(1:2, bar_vector, 1);
for dd = 1:length(pooled_tt_inf)
    plot(1:2, [pooled_tt_inf(dd) pooled_tt_sup(dd)], 'k');
    scatter(1, pooled_tt_inf(dd), 'filled', 'MarkerFaceColor', 'k');
    scatter(2, pooled_tt_sup(dd), 'filled', 'MarkerFaceColor', 'k');
end
scatter(1, ttinf_avg, 'filled', 'r');
scatter(2, ttsup_avg, 'filled', 'r');
errorbar(1:2, bar_vector, error_vector, '.');
text(1.5, max([ttinf_avg ttsup_avg]) + 0.025, sprintf('p = %0.3f', pooled_p));
xlabel('Condition')
xlim([0 3])
ylim([0 6])
ylabel('Mean total response time (s)')
xticks([1 2])
xticklabels({'Inf', 'Sup'})
title('Sequence task, pooled data')
axis square

%% Performance
perfmetric = 'Percent correct';     % metric for measuring performance: 'Percent correct' or 'LCS'
performance_inf = zeros(1, num_subjects);  % average performance in inferior field for each indiv
performance_sup = zeros(1, num_subjects);  % average performance in superior field for each indiv
difficulty_edges = [0 1.1 2.1 10];                % easy, medium, hard edges
diffcurve_inf = zeros(num_subjects, length(difficulty_edges)-1);         % performance vs difficulty for each field
diffcurve_sup = zeros(num_subjects, length(difficulty_edges)-1);     

for ii = 1:num_subjects
    expt = data(ii).expt;
    response = data(ii).response;
    
    sequence = expt.sequence(response.firsttime <= rt_thresh);
    field = expt.field(response.firsttime <= rt_thresh);
    response_seq = response.sequence(response.firsttime <= rt_thresh);
    
    % find proportion correct in every sequence
    prctcorrect = zeros(1, length(sequence));
    seq_difficulty = zeros(1, length(sequence));
    longest_match = zeros(1, length(sequence));
    for ss = 1:length(sequence)
        % find percent correct
        correct = sequence{ss} == response_seq{ss};
        prctcorrect(ss) = sum(correct)/length(correct);
        
        % find longest common sequence
        templatestr = sequencetostring(sequence{ss});
        responsestr = sequencetostring(response_seq{ss});
    
        [~, longest_match(ss), ~] = LCSubstr(templatestr, responsestr); 
        
        % calculate difficulty of the sequence
        u = length(unique(sequence{ss}));      % number of unique elements in the sequence
        c = max(diff([0 find(diff(sequence{ss})) numel(sequence{ss})]));       % max number of consecutive elements
        seq_difficulty(ss) = u/c;       
    end
    binned_difficulty = discretize(seq_difficulty, difficulty_edges);       % bin difficulties into 3 ranges (easy, medium, hard)
    
    switch perfmetric
        case 'Percent correct'
            performance = prctcorrect;
        case 'LCS'
            performance = longest_match;
    end
    
    % average performance and difficulty curves in each VF according to the desired performance metric 

    performance_inf(ii) = mean(performance(field == 0));
    performance_sup(ii) = mean(performance(field == 1));

    % difficulty curves
    for dd = 1:length(difficulty_edges)-1
        diffcurve_inf(ii, dd) = mean(performance(field == 0 & binned_difficulty == dd));   % inf field, curr difficulty
        diffcurve_sup(ii, dd) = mean(performance(field == 1 & binned_difficulty == dd));   % sup field, curr difficulty
    end
    
end

% pooled data
avg_performance_inf = mean(performance_inf);
se_performance_inf = std(performance_inf)/sqrt(num_subjects);
avg_performance_sup = mean(performance_sup);
se_performance_sup = std(performance_sup)/sqrt(num_subjects);

avg_diffcurve_inf = mean(diffcurve_inf, 1);
se_diffcurve_inf = std(diffcurve_inf, [], 1)/sqrt(num_subjects);
avg_diffcurve_sup = mean(diffcurve_sup, 1);
se_diffcurve_sup = std(diffcurve_sup, [], 1)/sqrt(num_subjects);

figure
hold on
bar([1 2], [avg_performance_inf avg_performance_sup]);
errorbar([1 2], [avg_performance_inf avg_performance_sup], [se_performance_inf se_performance_sup], 'LineStyle', 'none')
[~, p] = ttest(performance_inf, performance_sup);
xlim([0 3])
switch perfmetric
    case 'LCS'
        ylim([0 max([avg_performance_inf avg_performance_sup])+1]); 
        ylabel('Sequence report performance (LCS)')
    case 'Percent correct'
        ylim([0 1]);
        ylabel('Sequence report performance (% correct)')
end
title(sprintf('Average performance, p = %.2d', p))
xlabel('Visual field')
axis square

figure
hold on
errorbar(1:length(difficulty_edges)-1, avg_diffcurve_inf, se_diffcurve_inf, 'c')
errorbar(1:length(difficulty_edges)-1, avg_diffcurve_sup, se_diffcurve_sup, 'm')
for jj = 1:length(difficulty_edges)-1
    [~, p(jj)] = ttest(diffcurve_inf(:, jj), diffcurve_sup(:, jj));
end
xlim([0 length(difficulty_edges)])
switch perfmetric
    case 'LCS'
        ylim([0 max([avg_performance_inf avg_performance_sup])+1]); 
        ylabel('Sequence report performance (LCS)')
    case 'Percent correct'
        ylim([0 1]);
        ylabel('Sequence report performance (% correct)')
end
title('Average performance by difficulty');
xlabel('Sequence difficulty')


function sequencestring = sequencetostring(sequence)
    % input "sequence" = number sequence (array)
    % output "sequencestring" = input converted to equivalent string of letter directions
    
    sequencestring = '';
    
    for ss = 1:length(sequence)
        switch sequence(ss)
            case 0
                curr_ele = 'U';
            case 90
                curr_ele = 'L';
            case 180
                curr_ele = 'D';
            case 270
                curr_ele = 'R';
        end
        sequencestring = append(sequencestring, curr_ele);
    end
end
    


