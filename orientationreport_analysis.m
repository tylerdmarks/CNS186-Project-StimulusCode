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
offset = 0.2;
pooled_rt_inf = [];
pooled_rt_sup = [];
for ii = 1:num_subjects
    expt = data(ii).expt;
    response = data(ii).response;

    rt_inf = response.time(expt.field == 0);
    rt_sup = response.time(expt.field == 1);
    
    % exclude outliers
    rt_inf = rt_inf(rt_inf <= rt_thresh);
    rt_sup = rt_sup(rt_sup <= rt_thresh);
   
    figure(f0)
    violinplot(ii-offset, rt_inf', {sprintf('%d', ii)}, 'Width', 0.2, 'ViolinColor', [0.5 0.8 0.1], 'DataAlpha', 0.4, 'ViolinAlpha', 0.1);
    violinplot(ii+offset, rt_sup', {sprintf('%d', ii)}, 'Width', 0.2, 'ViolinColor', [0.3 0.2 0.9], 'DataAlpha', 0.4, 'ViolinAlpha', 0.1);
    p_rt(ii) = ranksum(rt_inf, rt_sup);
    text(ii, max([max(rt_inf) max(rt_sup)]), sprintf('p = %0.3f', p_rt(ii)));
    title('Response time')
    
    % pool average data across individuals
    pooled_rt_inf = [pooled_rt_inf mean(rt_inf)];
    pooled_rt_sup = [pooled_rt_sup mean(rt_sup)];
end
figure(f0)
xlim([0.5 num_subjects+0.5])
xlabel('Participant')
ylabel('Response time (s)')
% legend(plotdata, {'no primer', 'primer'})
title('Orientation task, by participant')
xticks(1:num_subjects)
for ss = 1:num_subjects
    labels{ss} = num2str(ss);
end
xticklabels(labels);

%%%% Plotting  pooled data


rtinf_avg = mean(pooled_rt_inf);
rtinf_se = std(pooled_rt_inf)/sqrt(length(pooled_rt_inf));
rtsup_avg = mean(pooled_rt_sup);
rtsup_se = std(pooled_rt_sup)/sqrt(length(pooled_rt_sup));

bar_vector = [rtinf_avg rtsup_avg];
error_vector = [rtinf_se rtsup_se];
[~, pooled_p] = ttest(pooled_rt_inf, pooled_rt_sup);

figure
hold on
% bar(1:2, bar_vector, 1);
for dd = 1:length(pooled_rt_inf)
    plot(1:2, [pooled_rt_inf(dd) pooled_rt_sup(dd)], 'k');
    scatter(1, pooled_rt_inf(dd), 'filled', 'MarkerFaceColor', 'k');
    scatter(2, pooled_rt_sup(dd), 'filled', 'MarkerFaceColor', 'k');
end
scatter(1, rtinf_avg, 'filled', 'r');
scatter(2, rtsup_avg, 'filled', 'r');
errorbar(1:2, bar_vector, error_vector, '.');
text(1.5, max([rtinf_avg rtsup_avg]) + 0.025, sprintf('p = %0.3f', pooled_p));
xlabel('Condition')
xlim([0 3])
ylim([0 3])
ylabel('Mean response time (s)')
xticks([1 2])
xticklabels({'Inf', 'Sup'})
title('Orientation task, pooled data')
axis square
%% Performance
performance_inf = zeros(1, num_subjects);  % average performance in inferior field for each indiv
performance_sup = zeros(1, num_subjects);  % average performance in superior field for each indiv
difficulty_edges = [0 0.26 0.49 1.1];                % easy, medium, hard edges
diffcurve_inf = zeros(num_subjects, length(difficulty_edges)-1);         % performance vs difficulty for each field
diffcurve_sup = zeros(num_subjects, length(difficulty_edges)-1);  
performance_inf_turn = zeros(1, num_subjects);
performance_inf_noturn = zeros(1, num_subjects);
performance_sup_turn = zeros(1, num_subjects);
performance_sup_noturn = zeros(1, num_subjects);

for ii = 1:num_subjects
    expt = data(ii).expt;
    response = data(ii).response;
    
    sequence = expt.sequence(response.time <= rt_thresh);
    field = expt.field(response.time <= rt_thresh);
    response_ori = response.orientation(response.time <= rt_thresh);
    
    % deterimine correct answer for each sequence
    % and determine sequence difficulty
    forward = [0 4 -4];       % results that equate to forward
    right = [1 5 -3 -7];        % results that equate to right
    behind = [2 6 -2 -6];
    left = [3 7 -1 -5];
    end_direction = zeros(1, length(sequence));         % for calculating end orientation
    seq_difficulty = zeros(1, length(sequence));        % sequence difficulty 
    turn = zeros(1, length(sequence));            % for tracking which sequences turn around the 180 degree point
    
    for ss = 1:length(sequence)
        curr_seq = sequence{ss};
        counter = 0;        % start a counter, 0 = facing forward
        for ee = 1:length(curr_seq)
            if curr_seq(ee) == 270
                counter = counter + 1;
            elseif curr_seq(ee) == 90
                counter = counter - 1;
            end
            
            if abs(counter) >= 3    
                turn(ss) = 1;
            end
        end
        if any(counter == forward)
            end_direction(ss) = 0;
        elseif any(counter == right)
            end_direction(ss) = 90;
        elseif any(counter == behind)
            end_direction(ss) = 180;
        elseif any(counter == left)
            end_direction(ss) = 270;
        end
        
        c = max(diff([0 find(diff(sequence{ss})) numel(sequence{ss})]));       % max number of consecutive elements
        seq_difficulty(ss) = 1/c;     % for now just directly scored by c
        
        
    end
    
    binned_difficulty = discretize(seq_difficulty, difficulty_edges);       % bin difficulties into 3 ranges (easy, medium, hard)

    correct_response = end_direction == response_ori;
    
    % average performance in each field
    performance_inf(ii) = sum(correct_response(field == 0))/sum(field == 0);
    performance_sup(ii) = sum(correct_response(field == 1))/sum(field == 1);
    
    % performance vs difficulty
    for dd = 1:length(difficulty_edges)-1
        diffcurve_inf(ii, dd) = sum(correct_response(field == 0 & binned_difficulty == dd))/sum(field == 0 & binned_difficulty == dd);   % inf field, curr difficulty
        diffcurve_sup(ii, dd) = sum(correct_response(field == 1 & binned_difficulty == dd))/sum(field == 1 & binned_difficulty == dd);   % sup field, curr difficulty
    end
    
    % performance in turnover trials vs no turnover
    performance_inf_turn(ii) = sum(correct_response(field == 0 & turn))/sum(field == 0 & turn);
    performance_inf_noturn(ii) = sum(correct_response(field == 0 & ~turn))/sum(field == 0 & ~turn);
    performance_sup_turn(ii) = sum(correct_response(field == 1 & turn))/sum(field == 1 & turn);
    performance_sup_noturn(ii) = sum(correct_response(field == 1 & ~turn))/sum(field == 1 & ~turn);
end

% pooled performance
avg_performance_inf = nanmean(performance_inf);
se_performance_inf = nanstd(performance_inf)/sqrt(num_subjects);
avg_performance_sup = nanmean(performance_sup);
se_performance_sup = nanstd(performance_sup)/sqrt(num_subjects);

% diffcurve
avg_diffcurve_inf = nanmean(diffcurve_inf, 1);
se_diffcurve_inf = nanstd(diffcurve_inf, [], 1)/sqrt(num_subjects);
avg_diffcurve_sup = nanmean(diffcurve_sup, 1);
se_diffcurve_sup = nanstd(diffcurve_sup, [], 1)/sqrt(num_subjects);

% as a function of turn
avg_performance_inf_turn = nanmean(performance_inf_turn);
se_performance_inf_turn = nanstd(performance_inf_turn)/sqrt(num_subjects);
avg_performance_sup_turn = nanmean(performance_sup_turn);
se_performance_sup_turn = nanstd(performance_sup_turn)/sqrt(num_subjects);
avg_performance_inf_noturn = nanmean(performance_inf_noturn);
se_performance_inf_noturn = nanstd(performance_inf_noturn)/sqrt(num_subjects);
avg_performance_sup_noturn = nanmean(performance_sup_noturn);
se_performance_sup_noturn = nanstd(performance_sup_noturn)/sqrt(num_subjects);


figure
hold on
bar([1 2], [avg_performance_inf avg_performance_sup]);
errorbar([1 2], [avg_performance_inf avg_performance_sup], [se_performance_inf se_performance_sup], 'LineStyle', 'none')
[~, p] = ttest(performance_inf, performance_sup);
xlim([0 3])
ylim([0 1])
title(sprintf('Average performance, p = %.2d', p))
xlabel('Visual field')
ylabel('Orientation report performance (% correct)')

figure
hold on
errorbar(1:length(difficulty_edges)-1, avg_diffcurve_inf, se_diffcurve_inf, 'c')
errorbar(1:length(difficulty_edges)-1, avg_diffcurve_sup, se_diffcurve_sup, 'm')
for jj = 1:length(difficulty_edges)-1
    [~, p(jj)] = ttest(diffcurve_inf(:, jj), diffcurve_sup(:, jj));
end
xlim([0 length(difficulty_edges)])
ylim([0 1])
title('Average performance by difficulty');
xlabel('Difficulty')
ylabel('Orientation report performance')

figure
hold on
bar([1 2], [avg_performance_inf_noturn avg_performance_sup_noturn], 1);
errorbar([1 2], [avg_performance_inf_noturn avg_performance_sup_noturn], [se_performance_inf_noturn se_performance_sup_noturn], 'LineStyle', 'none')
bar([3 4], [avg_performance_inf_turn avg_performance_sup_turn], 1);
errorbar([3 4], [avg_performance_inf_turn avg_performance_sup_turn], [se_performance_inf_turn se_performance_sup_turn], 'LineStyle', 'none')
[~, p_turn] = ttest(performance_inf_turn, performance_sup_turn);
[~, p_noturn] = ttest(performance_inf_noturn, performance_sup_noturn);
xlim([0 5])
ylim([0 1])
xlabel('Visual field')
ylabel('Orientation report performance (% correct)')

