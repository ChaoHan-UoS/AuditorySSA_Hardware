%%  Model Initialization
% The following code is adapted from that of Yarden's SSA model in auditory cortex
% This file generates a recurrent network and brings it to equilibrium

% close all
% clear all

all_plots = 0; % Controls all plots to display
plot_mixing = 0; % Plots the best frequencies for all neurons
plot_TC = 1; % Plots the shape of a typical tuning curve
plot_init = 0; % Plots the dynamics to reach steady state 

Conds{1} = 'Low';
Conds{2} = 'High';
Conds{3} = 'Equal';
Conds{4} = 'Diverse Broad';
Conds{5} = 'Deviant Alone F1';
Conds{6} = 'Deviant Alone F2';

nev_cond_code{1} = 'L';
nev_cond_code{2} = 'H';
nev_cond_code{3} = 'E';
nev_cond_code{4} = 'DB';
nev_cond_code{5} = 'DA1';
nev_cond_code{6} = 'DA2';

AXES_FONTSIZE = 10;
LineWidth = 1;
MarkerSize = 10;
COLORBAR_FONTSIZE = 10;


%% Parameters of the Model:
% Number of columns:
P = 5; %9

% Number of frequencies to which the network is sensitive
freq_stretch = 1; % Choose an integer >=1; this multiplies P to give M.
M = freq_stretch*P; % Thus, column best frequencies will be certain presentable frequencies.

ring_net = 0;

% Time constants:
tau = 0.001; %  membrane time constant of adaptive population (in seconds)
tau_E = 0.005; % excitatory neurons' membrane time constant (in seconds)
tau_I = 0.005; % inhibitory neurons' membrane time constant (in seconds) 
tau_a = 1; % adaptation time constant (in seconds)

% Connection strengths:
w_ee = 1.3; % exc2exc 1.4
w_ei = -1.2; % inh2exc -1.4 -1.4
w_ie = 0.75; % exc2inh 0.75 0.5
w_ii = -0.4; % inh2inh -0.5 0
% w_ee = 1.4; % exc2exc 1.4
% w_ei = -1.4; % inh2exc -1.4 -1.4
% w_ie = 0.5; % exc2inh 0.75 0.5
% w_ii = -0; % inh2inh -0.5 0
w_a = 0.5; % ad2exc
c = 20; % a_ss = c*Ea

w_ee1 = 0.075;  % exc2exc, between neighboring columns 0.045 0.170 0.080
w_ie1 = 0; % exc2inh, between neighboring columns
w_ee2 = 0; % exc2exc, from one column to its 2nd neighbor
w_ie2 = 0; % exc2inh, from one column to its 2nd neighbor

%Gain function of cortical cell:
Gain_thres = 0; % Threshold of gain 
Gain_slope = 2.5; % Slope of gain 2
Gain_max = 500;


%% Stimulus Parameters
t_eq = 5; % The time given to reach equilibrium (in seconds). It is important to allow enough time, otherwise the response to the first stimulus is distorted.
dt = 0.0001; % Time-step (in seconds)
post_stim = 0.080; % This is defined here so the simulations keeps runnng for 2*post_stim after the last stimulus offset

duration = 0.050; % Total duration of each stimulus (offset to onset) (in seconds) 0.050
ramp_dur = 0.005; % Durations of the ramps at the beginning and end of each stimulus (in seconds) 0.005
tmax = t_eq; % Maximum time to reach the equilibrium (in seconds) 
num_steps_eq = floor(tmax/dt); % Total number of steps to reach the equilibrium
        
%% Sensory Input
lambda_c = 2;
lambda_s_right = lambda_c;
lambda_s_left = lambda_c;

if ring_net
    Distances = toeplitz([0:floor(M/2) wrev(1:(floor(M-1)/2))]); % Distance of each frequency from all the others.
else
    Distances = toeplitz(0:(M-1));
end

Right_Dec_Args = Distances./lambda_s_right;
Left_Dec_Args = Distances./lambda_s_left;

curve_type = 'linear'; % Type of the tuning curve; see cases for options and details.

switch(curve_type)
    case 'exponential' % this one is the original type, used in Loebel et al. 2007. 
        h_outline = (triu(exp(-Right_Dec_Args)) + tril(exp(-Left_Dec_Args),-1)); % Magnitudes of input received by each column for each presented frequency
    
    case 'biased_exponential' % this is similar, but would be sharper and have inhibition from all farther off frequencies
        curve_bias = 1/3;
        h_outline = (1 + curve_bias)*(triu(exp(-Right_Dec_Args)) + tril(exp(-Left_Dec_Args),-1)) - curve_bias; % Ensures a magnitude of A at the tuning curve's peak
    
    case 'linear' % This type is triangular, following Eli's request that some columns do not receive input at all.
        h_outline = (triu(1-Right_Dec_Args) + tril(1-Left_Dec_Args,-1));
        h_outline(h_outline < 0) = 0;
        
    case 'mexican hat' % My own addition; has support in STRFs such as in deCharms et al. 1998; Loebel et al. 2007 state that lateral supression arises in the model due to synaptic depression 
        mex_dec = 1;
        mex_fre = 1;
        h_outline = (triu(cos(mex_fre*Right_Dec_Args).*exp(-mex_dec*Right_Dec_Args)) + tril(cos(mex_fre*Left_Dec_Args).*exp(-mex_dec*Left_Dec_Args),-1));

    case 'Gaussian' 
        h_outline = (triu(exp(-Right_Dec_Args.^2)) + tril(exp(-Left_Dec_Args.^2),-1)); % Magnitudes of input received by each column for each presented frequency
end

h_outline = h_outline(1:freq_stretch:end,:); % This gives the columns best frequencies that are spread evenly across the range.
h_outline(P+1:end,:) = []; % Cropping the rows of h_outline to the no. of columns, just in case.
h_outline = reshape(h_outline,[P 1 M]); % *) Aligning the tuning curves to the neurons.
h = h_outline;


%% Plotting the Tuning Curve 
if plot_TC && all_plots
    figure
    tc_column = 3; %floor(P/2);
    hold on
    plot(h_outline(tc_column,:), '-b.','LineWidth',LineWidth,'MarkerSize',MarkerSize)
    set(gca, 'FontSize', AXES_FONTSIZE,'YTick', [0 0.5 1],'XTick', [1 2 3 4 5])
    title(['Tuning Curve (Column ' num2str(tc_column) ')'])
    xlabel('Frequency')
    ylabel('Input as Fraction of the Input at the Best Frequency')
    xlim([1 5])
end

%% Plotting the Best Frequencies of All Neurons
[h_Peaks, Best_Freqs]= max(h,[],3);
if plot_mixing && all_plots
    figure
    imagesc(Best_Freqs)
    title('Best Frequencies of All Neurons')
    xlabel('No. of Neuron in Column')
    ylabel('Column')
    colorbar
end

%% Initializing Variables:
E = zeros(P,1); % Acitivity(firing rate) of all excitatory neurons (in Hz) 
I = zeros(P,1); % Activity of all inhibitory neurons (in Hz)
Ea = E; % Activity of neurons with adaptation (in Hz)
h_E = E;  % Sypnatic currents to excitatory population 
h_I = I; % Sypnatic currents to inhibitory population 
h_a = E; 
a = 0; % Amount of adaptation accumulated 

H_E = E; % % Sum of synaptic inputs to excitatory population at every time step
H_I = I; % % Sum of synaptic inputs to inhibitory population at every time step

E_act = zeros(P,floor(t_eq/dt)); % Activity of all excitatory neurons during the time allowed for reaching equilibrium
I_act = zeros(P,floor(t_eq/dt));
Ea_act = zeros(P,floor(t_eq/dt));
a_act = Ea_act;


%% Dynamic Loop: Initial run to obtain steady-state and find initial values
% The network is allowed to reach a steady-state with no sensory input. 

for i = 1:floor(t_eq/dt)
    % The dynamics of synaptic input to adaptative cell:
    h_a = h_a + (dt/tau)*(-h_a + 0);
    
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
    E(E >  Gain_max) = Gain_max;
    I = h_I;
    I = I - Gain_thres;
    I = I * Gain_slope;
    I(I <  0) = 0;
    I(I >  Gain_max) = Gain_max;

    % Tracking the activity of all neurons
    E_act(:,i) = E; 
    I_act(:,i) = I; 
    Ea_act(:,i) = Ea;
    a_act(:,i) = a;
    
end

%% Recording Equilibrium Conditions and Plot:
E_eq = E; 
h_E_eq = h_E;
I_eq = I; 
h_I_eq = h_I;
Ea_eq = Ea; 
h_a_eq = h_a;
a_eq = a;

if plot_init && all_plots
    figure
    plot(dt:dt:t_eq, E_act(3,:), 'b') % Firing rate of exc neuron in col 3
    hold on;
    plot(dt:dt:t_eq, I_act(3,:), 'r') % Firing rate of inh neuron in col 3
    plot(dt:dt:t_eq,Ea_act(3,:),'g') % Firing rate of adaptative neuron in col 3
    title('Mean Firing Rate of Column 3 during Initial Run')
    xlabel('time(s)')
    ylabel('Spikes/s')
    legend('E','I','Ea')
    
    figure
    plot(dt:dt:t_eq,a_act(3,:)) 
    title('adapation during Initial Run')
    xlabel('time(s)')
end

%% Saving
% E_act = [];% just to save memory
% I_act = [];
% Ea_act = [];
% a_act = [];
save('Simulation Results/initialization.mat');
disp('Initialization done');
