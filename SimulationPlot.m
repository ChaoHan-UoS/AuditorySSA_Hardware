%%  Plot simulation results 
close all
clear all


%% Propogation of PS
clear
load('Simulation Results/run_L_Tr1.mat');
% load('Simulation Results/run_L.mat');
t = Stim_Onsets(1);
time_win = 0.100; 

figure
imagesc([0 time_win*10^3], [1 P], E_act_overall(:,t:t+time_win/dt));
set(gca,'YTick',1:P,'FontSize',AXES_FONTSIZE)
xlabel('time(ms)')
ylabel('Column')
colormap jet
c = colorbar;
c.FontSize = COLORBAR_FONTSIZE;
c.Label.String = 'E (Spikes/s)';

%% Tuning curve 
clear
load('Simulation Results/run_L.mat')
figure
hold on
plot(h_outline(Rec_Column,:), '-k.','LineWidth',LineWidth,'MarkerSize',MarkerSize)
set(gca, 'FontSize', AXES_FONTSIZE,'YTick', [0 0.25 0.5 0.75 1], 'XTick', [1 2 3 4 5], 'XTickLabel', {'Q-2','Q-1','Q','Q+1','Q+2'})
title('Tuning Curve of Column Q')
xlabel('f')
ylabel('T_f^Q')
xlim([0 6])
y = [0.25 0; 0.25 0.25; 0 0; 0.25 0.75; 0.25 0];
b = bar(y,'FaceColor', 'flat', 'EdgeColor','flat');
legend(b,'Many-standards','Oddball')
legend('boxoff')
set(gca, 'FontSize', AXES_FONTSIZE);


%% Low condition 
clear
% load('Simulation Results/run_L.mat')
load(['Simulation Results/run_L_Tr1.mat']);
n_stim_plot = 8;
t = Stim_Onsets(1);
t_marg = 0.1; 
time_win = n_stim_plot*(duration+ISI)+t_marg; 
t_onset = t - floor(t_marg/dt); % Index of onset time 
t_offset = t + floor(time_win/dt) + 1; % Index of offset time 
tim = -t_marg:dt:time_win; 
t_step = 0.5;

% Oddball sequence
f = figure; 
% subplot(18,1,1);
% % newcolors = [0 1 0; 1 1 1; 1 0 0];
% % colororder(newcolors);
% % area(tim,[Spec_Temp(1,t_onset:t_offset,1)' 1-Spec_Temp(1,t_onset:t_offset,1)' Spec_Temp(1,t_onset:t_offset,2)'],'LineStyle','none');
% area(tim,Spec_Temp(1,t_onset:t_offset,1),'LineStyle','none','FaceColor',[0 0.4470 0.7410]);
% set(gca,'XTick',[0:t_step:time_win],'XTickLabel',[],'XLim',[-t_marg time_win],'FontSize', AXES_FONTSIZE)  
% axis([-t_marg time_win 0 4]);
% box off
% c = colorbar;
% c.Visible = 'off';

% subplot(18,1,2);
% area(tim,Spec_Temp(1,t_onset:t_offset,2),'LineStyle','none','FaceColor',[0.8500 0.3250 0.0980]);
% set(gca,'XTick',[0:t_step:time_win],'XTickLabel',[],'XLim',[-t_marg time_win],'FontSize', AXES_FONTSIZE)  
% axis([-t_marg time_win 0 4]);
% box off
% c = colorbar;
% c.Visible = 'off';

% set(gcf,'PaperUnits','normalized','PaperPosition',[0 0 1 1]) 
% set(gcf,'Units','normalized','position',get(gcf,'PaperPosition')) 
% saveas(f,'figure/Odd_Seq.pdf');


% Adaptive population
subplot(7,1,5); % Dynamics of adaptation variable a in 3 cols
plot(tim,a_act_overall(Rec_Column-Freq_Diff/2, t_onset:t_offset),'Color',[0 0.4470 0.7410],'LineWidth',LineWidth,'MarkerSize',MarkerSize); % Std column, blue
hold on;
plot(tim,a_act_overall(Rec_Column, t_onset:t_offset),'Color',[0.4940 0.1840 0.5560],'LineWidth',LineWidth,'MarkerSize',MarkerSize); % Middle column, purple
plot(tim,a_act_overall(Rec_Column+Freq_Diff/2, t_onset:t_offset),'Color',[0.8500 0.3250 0.0980],'LineWidth',LineWidth,'MarkerSize',MarkerSize); % Dev column, red
set(gca,'XTick',[0:t_step:time_win],'XTickLabel',[],'XLim',[-t_marg time_win],'YLim',[0 15],'TickDir','out','FontSize', AXES_FONTSIZE)  
ylabel('Adaptive Thres (Hz)')
box off
c = colorbar;
c.Visible = 'off';

subplot(7,1,6); % Dynamics of adaptive population activities A_a in 3 cols
plot(tim,Ea_act_overall(Rec_Column-Freq_Diff/2, t_onset:t_offset),'Color',[0 0.4470 0.7410],'LineWidth',LineWidth,'MarkerSize',MarkerSize);
hold on;
plot(tim,Ea_act_overall(Rec_Column, t_onset:t_offset),'Color',[0.4940 0.1840 0.5560],'LineWidth',LineWidth,'MarkerSize',MarkerSize);
plot(tim,Ea_act_overall(Rec_Column+Freq_Diff/2, t_onset:t_offset),'Color',[0.8500 0.3250 0.0980],'LineWidth',LineWidth,'MarkerSize',MarkerSize);
set(gca,'XTick',[0:t_step:time_win],'XTickLabel',[],'XLim',[-t_marg time_win],'YLim',[0 15],'TickDir','out','FontSize', AXES_FONTSIZE)  
ylabel('A_a (Hz)')
box off
legend('Std','Rec','Dev','Box','off')
c = colorbar;
c.Visible = 'off';

subplot(7,1,7); % Propogation of adaptive population activities A_a across 5 cols
image([-t_marg time_win], [1 5], Ea_act_overall(:,t_onset:t_offset));
set(gca, 'XTick',[0:t_step:time_win],'XLim',[-t_marg time_win],'YTick',1:5,'TickDir','out','FontSize', AXES_FONTSIZE);
ylabel('Column');
box off
colormap(gca, viridis(15));
c = colorbar;
c.TickDirection = 'out';
c.TickLength  = 0.03;
c.Box = 'off';
c.FontSize = COLORBAR_FONTSIZE;
c.Label.String = 'A_a (Hz)';
c.Limits = [0 15];
c.Ticks = [0 5 10 15];

set(gcf,'PaperUnits','normalized','PaperPosition',[0 0 1 1]) 
set(gcf,'Units','normalized','position',get(gcf,'PaperPosition')) 
saveas(f,'figure/Adap_popu.pdf');

% E/I population
f = figure;
subplot(7,1,1);  % Dynamics of E/I population activities A_e/A_i in the mid col
plot(tim,E_act_overall(Rec_Column, t_onset:t_offset),'Color','#EA6B66',...
     'LineWidth',LineWidth,'MarkerSize',MarkerSize);
hold on;
plot(tim,I_act_overall(Rec_Column,t_onset:t_offset),'Color','#7EA6E0',...
     'LineWidth',LineWidth,'MarkerSize',MarkerSize);
set(gca,'XTick',[0:t_step:time_win],'XTickLabel',[],'XLim',[-t_marg time_win],'TickDir','out','FontSize', AXES_FONTSIZE)  
ylabel('Hz')
box off
legend('Exc','Inh','Box','off')
c = colorbar;
c.Visible = 'off';

subplot(7,1,2); % Propogation of Exicitatory population activities A_e across 5 cols
image([-t_marg time_win], [1 5], E_act_overall(:,t_onset:t_offset));
set(gca, 'XTick',[0:t_step:time_win],'XLim',[-t_marg time_win],'YTick',1:5,'TickDir','out','FontSize', AXES_FONTSIZE);
ylabel('Column');
% xlabel('time(s)')
box off
colormap(gca, viridis(200));
c = colorbar;
c.TickDirection = 'out';
c.TickLength  = 0.03;
c.Box = 'off';
c.FontSize = COLORBAR_FONTSIZE;
c.Label.String = 'A_e (Hz)';
c.Limits = [0 200];
c.Ticks = [0 50 100 150 200];

set(gcf,'PaperUnits','normalized','PaperPosition',[0 0 1 1]) 
set(gcf,'Units','normalized','position',get(gcf,'PaperPosition')) 
saveas(f,'Figure/EI_Popu.pdf');



%% Many-standards
clear
load('Simulation Results/run_DB.mat')
n_stim_plot = 8;
t_onset = Stim_Onsets(1);
t_offset = t_onset + ceil(n_stim_plot*(duration+ISI)/dt) - 1;
tim = dt:dt:n_stim_plot*(duration+ISI);

figure % Adaptation
subplot(4,1,1);
hold on
% plot(tim,Ea_act_overall(Rec_Column - 2*Freq_Diff/2, t_onset:t_offset),'c','LineWidth',LineWidth);
% plot(tim,Ea_act_overall(Rec_Column - Freq_Diff/2, t_onset:t_offset),'b','LineWidth',LineWidth);
plot(tim,Ea_act_overall(Rec_Column, t_onset:t_offset),'g','LineWidth',LineWidth);
plot(tim,Ea_act_overall(Rec_Column + Freq_Diff/2, t_onset:t_offset),'r','LineWidth',LineWidth);
plot(tim,Ea_act_overall(Rec_Column + 2*Freq_Diff/2, t_onset:t_offset),'m','LineWidth',LineWidth);
xlim([0 n_stim_plot*(duration+ISI)])
set(gca, 'XTickLabel',[],'FontSize',AXES_FONTSIZE);
ylabel('Spikes/s')
legend('Std2','Std1','Rec','Dev1','Dev2')
c = colorbar;
c.Visible = 'off';
title('Adaptation Layer')

subplot(4,1,2);
hold on
% plot(tim,a_act_overall(Rec_Column - 2*Freq_Diff/2, t_onset:t_offset),'c','LineWidth',LineWidth);
% plot(tim,a_act_overall(Rec_Column - Freq_Diff/2, t_onset:t_offset),'b','LineWidth',LineWidth);
plot(tim,a_act_overall(Rec_Column, t_onset:t_offset),'g','LineWidth',LineWidth);
plot(tim,a_act_overall(Rec_Column + Freq_Diff/2, t_onset:t_offset),'r','LineWidth',LineWidth);
plot(tim,a_act_overall(Rec_Column + 2*Freq_Diff/2, t_onset:t_offset),'m','LineWidth',LineWidth);
xlim([0 n_stim_plot*(duration+ISI)])
set(gca, 'XTickLabel',[],'FontSize',AXES_FONTSIZE);
ylabel('a')
c = colorbar;
c.Visible = 'off';

subplot(4,1,3);
imagesc([0 n_stim_plot*(duration+ISI)], [1 4], reshape(Spec_Temp(1,t_onset:t_offset,:),length(t_onset:t_offset),4)');
map = [1 1 1; 0 0 0];
colormap(gca, map)
set(gca, 'YTick', [1 2 3 4], 'XTickLabel', [], 'YTickLabel', [1 2 4 5]);
c = colorbar;
c.Visible = 'off';

subplot(4,1,4);
imagesc([0 n_stim_plot*(duration+ISI)], [1 5], Ea_act_overall(:,t_onset:t_offset));
set(gca, 'YDir', 'normal', 'YTick', [1 2 3 4 5], 'FontSize', AXES_FONTSIZE);
ylabel('Column');
colormap(gca, hot);
c = colorbar;
c.FontSize = COLORBAR_FONTSIZE;
c.Label.String = 'E (Spikes/s)';

figure % E/I
subplot(3,1,1);
plot(tim,E_act_overall(Rec_Column,t_onset:t_offset),'b','LineWidth',LineWidth);
hold on
xlim([0 n_stim_plot*(duration+ISI)])
set(gca, 'XTickLabel',[],'FontSize',AXES_FONTSIZE);
ylabel('E(Spikes/s)')
c = colorbar;
c.Visible = 'off';
title('Cortex')

subplot(3,1,2);
imagesc([0 n_stim_plot*(duration+ISI)], [1 4], reshape(Spec_Temp(1,t_onset:t_offset,:),length(t_onset:t_offset),4)');
map = [1 1 1; 0 0 0];
colormap(gca, map)
set(gca,'YDir', 'normal', 'YTick', [1 2 3 4], 'XTickLabel', [], 'YTickLabel', [1 2 4 5]);
c = colorbar;
c.Visible = 'off';

subplot(3,1,3);
imagesc([0 n_stim_plot*(duration+ISI)], [1 5], E_act_overall(:,t_onset:t_offset));
set(gca, 'YDir', 'normal', 'YTick', [1 2 3 4 5], 'FontSize', AXES_FONTSIZE);
ylabel('Column');
colormap(gca, hot);
c = colorbar;
c.FontSize = COLORBAR_FONTSIZE;
c.Label.String = 'E (Spikes/s)';


%% SSA in other protocols
clear
load('Simulation Results/meta_data.mat'); 
Conds_Arr_plot = [2 1 3 6 4]; % Odd Standard in Col4, Odd deviant in Col4, Equal, Dev alone in Col4, Dev among std in Col4 [2 1 3 6 4]
time_win = 0.100;  % (in seconds) Spiking response within a window to plot
E_mean = zeros(length(Conds_Arr), ceil(time_win/dt));
Sp_m_arr = zeros(length(Conds_Arr),1); % mean of spike count of all protocols
Sp_sd_arr = Sp_m_arr; % standard deviation of spike count of all protocols
Sp_se_arr = Sp_m_arr; % standard error of the mean of all protocols
jj = 1; % Counter of number of protocols

% Compute dev alone first for normalization
j = 0; % Number of deviant reponses in column 4 
Spcount = zeros(1, n_stim*num_trials);
for trnum = 1:num_trials
    load(['Simulation Results/run_' nev_cond_code{6} '_Tr' num2str(trnum) '.mat']);
    for ns = 1:n_stim
        if Oddball(ns) == Rec_Column+1
            Spcount(j+1) = sum(E_act_overall(Rec_Column, Stim_Onsets(ns):Stim_Onsets(ns)+ceil(time_win/dt)-1)*dt);
            j = j + 1;
        end
    end
end
Spcount = Spcount(1,1:j); % Only deviant responses
Sp_m_DevA = mean(Spcount);

for ii = Conds_Arr_plot % Choose specific protocols
    j = 0; % Number of stimuli to column 4 in each protocol
    Spcount = zeros(1, n_stim*num_trials);
    for trnum = 1:num_trials
        load(['Simulation Results/run_' nev_cond_code{ii} '_Tr' num2str(trnum) '.mat']);
%         load(['Simulation Results/run_' nev_cond_code{ii} '.mat']);
        for ns = 1:n_stim
            if Oddball(ns) == Rec_Column+1
                E_mean(jj,:) = E_mean(jj,:) + E_act_overall(Rec_Column, Stim_Onsets(ns):Stim_Onsets(ns)+ceil(time_win/dt)-1);
                Spcount(j+1) = sum(E_act_overall(Rec_Column, Stim_Onsets(ns):Stim_Onsets(ns)+ceil(time_win/dt)-1)*dt);
                j = j + 1;
            end
        end 
    end
    E_mean(jj,:) = E_mean(jj,:)/j;
    Spcount = Spcount(1,1:j); % Only deviant responses
    Spcount = Spcount/Sp_m_DevA;
    Sp_m = mean(Spcount);
    Sp_sd = std(Spcount);
    Sp_se = Sp_sd/sqrt(n_stim*Probs);
    Sp_m_arr(jj) = Sp_m; % Mean of spike count for bar graph
    Sp_sd_arr(jj) = Sp_sd; % Standard deviation of spike count
    Sp_se_arr(jj) = Sp_se; % Standard error of spike count
    jj = jj + 1;
end

% SI and CSI computation
SI = (Sp_m_arr(2) - Sp_m_arr(1))/(Sp_m_arr(2) + Sp_m_arr(1))
CSI = (Sp_m_arr(2) - Sp_m_arr(5))/(Sp_m_arr(2) + Sp_m_arr(5))

% Plotting
color{1} = [0 0.4470 0.7410]; % Blue
color{2} = [0.8500 0.3250 0.0980]; % Red
color{3} = [0.3010 0.7450 0.9330]; % Cyan
color{4} = [0.4660 0.6740 0.1880]; % Green
color{5} = [0.9290 0.6940 0.1250]; % Yellow

f = figure; % Mean firing rate of each protocol
subplot(2,2,1);
tim = dt:dt:time_win; % (in seconds)
tim = tim*1000; % (in milliseconds)
hold on
for  i = 1:length(Conds_Arr)
    plot(tim, E_mean(i,:), 'Color',color{i},'LineWidth',LineWidth)
end
set(gca, 'FontSize', AXES_FONTSIZE,'XTick', [0 50 100]);
xlabel('time(ms)')
ylabel('A_e(Spikes/s)')
legend('Standard in OD','Deviant in OD','Equal','Dev alone','Dev among std')
legend('boxoff')
axis square

set(gcf,'PaperUnits','normalized','PaperPosition',[0 0 1 1]) 
set(gcf,'Units','normalized','position',get(gcf,'PaperPosition')) 
saveas(f,'Figure/Mean_f.pdf')

f = figure;
subplot(3,3,1);
b = bar(Sp_m_arr, 'FaceColor', 'flat', 'EdgeColor','flat');
for n = 1:length(Conds_Arr)
    b.CData(n,:) = color{n};
end
axis equal
ylabel('Spike Count');
xticklabels({'Standard in OD','Deviant in OD','Equal','Dev alone','Dev among std',})
xtickangle(45)
set(gca,'XLim',[0 6],'YLim',[0 1.1],'FontSize', AXES_FONTSIZE);
box off
hold on
er = errorbar(Sp_m_arr,Sp_se_arr); 
er.LineWidth = LineWidth;     
er.Color = [0 0 0];
er.LineStyle = 'none';  
axis square

set(gcf,'PaperUnits','normalized','PaperPosition',[0 0 1 1]) 
set(gcf,'Units','normalized','position',get(gcf,'PaperPosition')) 
saveas(f,'Figure/Bar.pdf');
disp('done')



% figure
% barplot = zeros(1,6);
% [barplot I] = sort(barplot);
% barplot = barplot/max(barplot);
% b = bar(barplot);
% b.FaceColor = 'flat';
% for n = 1:6
%     b.CData(n,:) = colorb{I(n)};
% end

