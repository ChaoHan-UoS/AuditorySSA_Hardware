%% Model Run
% This code runs the network saved and brought to equilibrium by model_init, according to the specified nev_cond

% nev_cond = 1;
% A = 1;
% ISI = 1;
% Freq_Diff = 1;
% Probs_Arr = 1;
% Probs_F2L = 1;
% Rec_Column = 1;

tic
t_prot = n_stim*(duration + ISI) + 2*post_stim; % % Total time of the protocol
tmax_tot = t_eq + t_prot; % Maximum time that the simulation will reach (in seconds) 
num_steps = floor(tmax_tot/dt); % Total number of steps in the simulation

Freq_Scale = 1:M;

Prob.L = Probs_F2L;
Freqs_Pres.L = [Rec_Column - Freq_Diff/2; Rec_Column + Freq_Diff/2]; %Rec_Column = 3; Freq_Diff  = 2; [2; 4]
    
Prob.H = 1 - Probs_F2L;
Freqs_Pres.H = [Rec_Column + Freq_Diff/2; Rec_Column - Freq_Diff/2]; %[4; 2]

Prob.E = 0.5;
Freqs_Pres.E = [Rec_Column - Freq_Diff/2; Rec_Column + Freq_Diff/2]; %[2; 4]
        
Prob.DB = Probs_F2L;
if Prob.DB == 0.25
%     flank_freq = (1/Probs_F2L - 2)/2;
%     Freqs_Pres.DB = sort([(Rec_Column - Freq_Diff/2) - (0:Freq_Diff:flank_freq*Freq_Diff), (Rec_Column + Freq_Diff/2) + (0:Freq_Diff:flank_freq*Freq_Diff)]);
%     Freqs_Pres.DB = Freqs_Pres.DB';
    %[2 4 6 8 10 12 14 16 18 20]'
    Freqs_Pres.DB = [1 2 4 5];
%     Freqs_Pres.DB = [4 5 7 8];
%     Freqs_Pres.DB = [2 4  6 8];
%     Freqs_Pres.DB = [3 4  6 7];
    Freqs_Pres.DB = Freqs_Pres.DB';
end   
    
Prob.DA1 = Probs_F2L;
Freqs_Pres.DA1 = Rec_Column - Freq_Diff/2; %2
        
Prob.DA2 = Probs_F2L;
Freqs_Pres.DA2 = Rec_Column + Freq_Diff/2; %4

switch nev_cond
    case('Low')
        Probs = Prob.L;
        Freqs_Pres = Freqs_Pres.L;
    case('High')
        Probs = Prob.H;
        Freqs_Pres = Freqs_Pres.H;
    case('Equal')
        Probs = Prob.E;
        Freqs_Pres = Freqs_Pres.E;
    case('Diverse Broad')
        Probs = Prob.DB;
        Freqs_Pres = Freqs_Pres.DB;
    case('Deviant Alone F1')
        Probs = Prob.DA1;
        Freqs_Pres = Freqs_Pres.DA1;
    case('Deviant Alone F2')
        Probs = Prob.DA2;
        Freqs_Pres = Freqs_Pres.DA2;   
end

% Indices of presented frequencies in frequency scale
Freqs_pres_inds = zeros(1,length(Freqs_Pres)); 
for frin = 1:length(Freqs_Pres)
    Freqs_pres_inds(frin) = find(Freq_Scale == Freqs_Pres(frin)); 
end

h_sq = h(:,:,Freqs_pres_inds); % Value of turing curve of each neuron at frequencies in certain protocal 

%% Oddball Sequence Generation
switch nev_cond
    case('Low')
        rng(666);
        Oddball = [Freqs_Pres(2)*ones(1,ceil(Probs*n_stim)) Freqs_Pres(1)*ones(1,ceil((1-Probs)*n_stim))]; % Column 4 is deivant and Column 2 standard
        Oddball = Oddball(randperm(n_stim));
%         Oddball = [2,2,2,4,2,2,4,2,2,2,4,2,2,2,4,2,2,4,2,4,2,2,4,2,2,2,2,2,4,2,4,2,2,2,2,2,2,2,4,2];
    case('High')
        rng(666);
        Oddball = [Freqs_Pres(1)*ones(1,ceil(Probs*n_stim)) Freqs_Pres(2)*ones(1,ceil((1-Probs)*n_stim))];
        Oddball = Oddball(randperm(n_stim));
    case('Equal')
        rng(666);
        Oddball = [Freqs_Pres(2)*ones(1,ceil(Probs*n_stim)) Freqs_Pres(1)*ones(1,ceil((1-Probs)*n_stim))];
        Oddball = Oddball(randperm(n_stim));
    case('Diverse Broad')
        rng(666);
        freq_reps = n_stim/length(Freqs_Pres);
        Oddball = repmat(Freqs_Pres,1,freq_reps);
        Oddball = Oddball(randperm(n_stim));
    case('Deviant Alone F1')
        rng(666);
        Oddball = [Freqs_Pres(1)*ones(1,ceil(Probs*n_stim)) zeros(1,ceil((1-Probs)*n_stim))];
        Oddball = Oddball(randperm(n_stim));
    case('Deviant Alone F2')
        rng(666);
        Oddball = [Freqs_Pres(1)*ones(1,ceil(Probs*n_stim)) zeros(1,ceil((1-Probs)*n_stim))];
        Oddball = Oddball(randperm(n_stim));
end

Spec_Temp = zeros(1,num_steps,length(Freqs_pres_inds)); % Spectro-temporal structure of the oddball stimulus sequence Xi_f(t)
ramp_step = dt/ramp_dur;
Single_Stim = [ramp_step:ramp_step:1, ones(1,floor((duration - 2*ramp_dur)/dt)), wrev(ramp_step:ramp_step:1)]; % Waveform of each single stimulus
Stim_Onsets = zeros(1,n_stim); % Onset time indices of each stimulus (no.1 ... no.100)
Time_Ind = zeros(n_stim,2); % Onset & offset time indices of each stimulus

for ns = 1:n_stim    
    Time_Ind(ns,:) = floor(t_eq/dt) + floor((ISI+duration)/dt)*(ns - 1) + [0, (length(Single_Stim)-1)]; 
    if Oddball(ns)
        Spec_Temp(1, Time_Ind(ns,1):Time_Ind(ns,2), Freqs_Pres == Oddball(ns)) = Single_Stim; % Amplitude of each stimulus along time indices
    end
    Stim_Onsets(ns) = Time_Ind(ns,1);
end

%% Returning to Equilibrium Conditions:
E = E_eq;
I = I_eq;
a = a_eq; 
h_a = h_a_eq;

% The following lines are for tracking activity of all neurons (made into comments for memory limitations)
E_act_overall = zeros(P,num_steps); 
E_act_overall(:,1:floor(t_eq/dt)) = E_act(:,:);
I_act_overall = zeros(P,num_steps); 
I_act_overall(:,1:floor(t_eq/dt)) = I_act(:,:);
Ea_act_overall = zeros(P,num_steps); 
Ea_act_overall(:,1:floor(t_eq/dt)) = Ea_act(:,:);
a_act_overall = zeros(P,num_steps); 
a_act_overall(:,1:floor(t_eq/dt)) = a_act(:,:);

%% Dynamic Loop

i = floor(t_eq/dt) + 1;
while i < num_steps
    % Adding uncertainty to network paramters
    tau_a = tau_a_op + tau_a_op*p_fluc*(rand-0.5)*2;
    w_ee = w_ee_op + w_ee_op*p_fluc*(rand-0.5)*2;
    w_ei = w_ei_op + w_ei_op*p_fluc*(rand-0.5)*2;
    w_ie = w_ie_op + w_ie_op*p_fluc*(rand-0.5)*2;
    w_ii = w_ii_op + w_ii_op*p_fluc*(rand-0.5)*2;
    w_a = w_a_op + w_a_op*p_fluc*(rand-0.5)*2;
    c = c_op + c_op*p_fluc*(rand-0.5)*2;
    w_ee1 = w_ee1_op + w_ee1_op*p_fluc*(rand-0.5)*2;
    Gain_slope = Gain_slope_op + Gain_slope_op*p_fluc*(rand-0.5)*2;
    
    % Calculating the sensory input:
    s_E = A*repmat(Spec_Temp(1,i,:),P,1).*h_sq;
    
    % The dynamics of synaptic input to adaptative cell:
    h_a = h_a + (dt/tau)*(-h_a + sum(s_E,3));
    
    % Implementing the non-linearity:
    Ea = h_a - a;
    Ea(Ea <  0) = 0;
    
    % The dynamics of adaptation variable:
    a = a + (dt/tau_a)*(-a + c*Ea);
    
    % Pre-calculation for the inter-column exc2exc gain:
    E1 = [E(2:P); ring_net*E(1)] + [ring_net*E(P); E(1:P-1)]; % Excitatory input from nearest neighboring column
    E2 = [E(3:P); ring_net*E(1:2)] + [ring_net*E(P-1:P); E(1:P-2)]; % Excitatory input from 2nd-nearest neighboring column
    
    % The total synaptic input cell receives:
    H_E = w_ee*E + w_ei*I + w_ee1.*E1 + w_ee2*E2 + w_a*Ea;
    H_I = w_ie*E + w_ii*I + w_ie1*E1 + w_ie2.*E2;
    
    % The dynamics of synaptic input to exc and inh cells:
    h_E = h_E + (dt/tau_E)*(-h_E + H_E);
    h_I = h_I + (dt/tau_I)*(-h_I + H_I);
    
    % Implementing the non-linearity:
    E = h_E;
    E = E - Gain_thres;
    E = E * Gain_slope;
    E(E <  0) = 0;
    I = h_I;
    I = I - Gain_thres;
    I = I * Gain_slope;
    I(I <  0) = 0;
    
    % Tracking activity:
    E_act_overall(:,i) = E;
    I_act_overall(:,i) = I;
    Ea_act_overall(:,i) = Ea;
    a_act_overall(:,i) = a;
    
    i = i + 1;
end

%% Saving
if save_results
    filename2 = ['Simulation Results/run_' nev_cond_code{nev_co}];
    
    if length(Probs_Arr) > 1
        filename2 = [filename2 '_P' num2str(Probs_F2L)];
    end
    if length(Freq_Diff_Arr) > 1
        filename2 = [filename2 '_F' num2str(Freq_Diff)];
    end
    if length(ISI_Arr) > 1
        filename2 = [filename2 '_ISI' num2str(ISI)];
    end
    if length(A_Arr) > 1
        filename2 = [filename2 '_A' num2str(A)];
    end
    
    filename2 = [filename2 '_Tr' num2str(tr)];
    filename2 = [filename2 '.mat'];
    save(filename2)
    
end

curr_time2 = clock;
disp([num2str(curr_time2(4)) ':' num2str(curr_time2(5)) ' ' nev_cond ' condition, ISI = ' num2str(ISI) 's, A = ' num2str(A)  'Hz, Trial ' num2str(tr) ' done']) 

toc
