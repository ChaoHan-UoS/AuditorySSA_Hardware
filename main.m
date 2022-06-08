%% Run 6 conditions
% The code is adapted from Yarden TS, Nelken I (2017)

close all
clear

Init = 1; 
if Init
    run model_init.m
else
    load(['Simulation Results/initialization.mat'])
end

Rec_Column = 3; %5
Freq_Diff_Arr = [2];

% Protocol Parameters:
Conds_Arr = [1 2 3 4 6]; % 6 basic conditions;  [1 2 3 4 6]
Probs_Arr = [0.25]; % Sets the probabilitis of F2 in the Low Condition
A_Arr = [15]; % Amplitude of stimuli
ISI_Arr = [0.300]; % Interval between stimuli (offset to onset)

n_stim = 800; % Total no. of stimuli in each trial (Best take a product of 4), 100
num_trials = 1;

for nevco = 1:length(Conds_Arr) % There are 6 basic conditions;
    for prob = 1:length(Probs_Arr)
        for tr = 1:num_trials
            for cc = 1:length(Freq_Diff_Arr)
                for bb = 1:length(ISI_Arr)
                    for aa = 1:length(A_Arr)
                        
                        save('Simulation Results/meta_data.mat', 'Rec_Column', 'Freq_Diff_Arr', ...
                        'Conds_Arr', 'Probs_Arr', 'A_Arr', 'ISI_Arr', 'n_stim', 'num_trials', 'nevco', 'prob', 'tr', 'cc', 'bb', 'aa', 'dt', 'nev_cond_code')
                        clear
                        load('Simulation Results/meta_data.mat')
%                         run model_init.m
                        load('Simulation Results/initialization.mat')
                        
                        nev_co  = Conds_Arr(nevco);
                        nev_cond  = Conds{nev_co};
                        Probs_F2L = Probs_Arr(prob);
                        Freq_Diff = Freq_Diff_Arr(cc);
                        A = A_Arr(aa);
                        ISI = ISI_Arr(bb);
                        
                        save_results = 1;
                        run model_run.m
                    end
                end
            end
        end
    end
    disp(['Calculation for ' Conds{nev_co} ' condition completes'])
    disp('--------------------------------------------------------------------------------------')
end


