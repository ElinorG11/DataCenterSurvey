% The purpose of this code is to simulate frequency behavior during LLM 
% training. We're particularly interested in how the transient power peak 
% at the initiation of training (represented as a step function) might  
% cause instability in the frequency response. Although the training process 
% also involves interesting events like GPU communication and checkpointing, 
% this particular simulation focuses on the effects of the initial power peak.

% The dynamic equations we aim to simulate are:
% (d/dt) delta = Δw
% (d/dt) Δw = K(Pref - PL(t)) - (3K|Eg||E|)/X * sin(delta) - (K/D) * Δw

% We found the equilibrium point
% w* = ws
% delta* = sin^-1(X/(3|Eg||E|) * (Pref-PL*)
% Pg* = PL* - Pref
% P* = Pref

% Now we wish to find the scale of the parameters:
% First, it is clear that at t=0 we are at steady state and PL* = 0 -> 
% Pref = 0.
% Next, examine the expression X/(3|Eg||E|). 
% E is the generator's EMF in absolute value in units of [V]
% Eg is the infinite bus voltage in absolute value in units of [V]
% X is the Reactance of the transmission line in units of [Ohm]

% The units are [Ohm] / [V^] = 1 / [W]
% Thus, we denote 1/Px = X/(3|Eg||E|). Note that Px represents the power
% that may be passed therough the distribution network.
% Now, delta* = sin^-1((1/Px) * (0 - PL*)).
% So, for a smart planning of the system, we would choose the scale of the
% productino from the generator, in the same scale as the power that may be
% transmitted through the distribution network. i.e., about Px approx 2PL*
% Following, we would like to choose K/D. We know that D/K is associated
% with the time decay of the system, and has units of [s]. Thus, K/D, as
% seen from the calculation in the file "" is about 1 [s]. This calculation
% is based on the inherent parametersassociated with the mechanical 
% properties of generators.
% Finally, we want to determine the scale of K. The units are [1/W*s^2],
% thus we will choose it proportionally to 1/Px. Under these assumptions we
% will result in logical values for the system we are simulating.

% Using the definitions above, the dynamic equations become:
% (d/dt) delta = Δw,
% (d/dt) Δw = K *(Pref - PL(t)) - K * Px * sin(delta) - alpha * Δw,
% where alpha = K/D.

clear; clc;

% Simulation parameters
PLstep = 50e6; % This is the step value (the maximal power drawn by the data 
               % center) in units of [W]
Px = 2*PLstep; % This parameter represents the expression (3|Eg||E|)/x 
               % where X, |Eg| and |E| are defined above.
Prt = Px; % This is the rated power of the generator. A sensible value is 
          % about twice the size of the demand.
fs = 60; % Nominal electrical frequency in units of [Hz]
ws = 2 * pi * fs; % Nominal electrical frequency in units of [rad/s]
K = 2.2e-04 * ws^2/Prt; % Generator's inertia constant in units of [1/(W*s^2)];
                   % Note that the numeric value of 30 has no units.
Pref = 0.5*PLstep; % 3 phase generator's reference power in units of [W] 
                   % Note: if we set Pref to be positive, i.e., at the
                   % initiation of the process, when PL=0, the generator
                   % will inject power into the network (since Pref>0).
                   % Thus, the angle delta will be poisitive as well (delta
                   % >0).
alpha = 100; % This parameter represents the fraction K/D in units of [1/s]. 
             % D is the generator's droop constant in units of [1/(W*s)]. 
             % Note that by uncreasing the value of alpha -> we decrease D
             % -> We increase the damping (because in the dynamic equations
             % of the generator it is 1/D).
% delta is the electrical angle of the generator in units of [rad].
% P is the 3 phase generator's power output in units of [W].

% Simulink parameters
SimTime = 10; % Simulation time in units of [s]
RelTol = 1e-4; % Simulation accuracy. Based on trial & error.
MaxStep = 1e-3; % Simulation max step size in units of [s] 
                % Based on trial & error.

% Simulink
disp('Running Simulink.');
disp('Please wait...');
DataCenterSim; % Open simulink
sim(bdroot); % Run Simulink
disp('Done running Simulink.');

% Plot results
% In these figures we see two interesting phenomena. The first, is the
% phenomena we excpected. We see the instability in the frequency around
% the time of the step in the power of the load. The second, is by
% experimenting with different values of D (basically by altering alpha).
% We see that by increasing the damping (decreasing D and increasing alpha)
% we can damp the oscillations in the frequency and in the angle delta 
% much faster, but on account of the required power, since we will need to
% supply more power for longer time.
figure(1);
subplot(4,1,1);
plot(ts, PL/1e6);
hold on;
ylabel('PL [MW]');

subplot(4,1,2);
plot(ts, P/1e6);
hold on;
ylabel('P [MW]');

subplot(4,1,3);
plot(ts, (DeltaOmega+ws)/(2*pi));
hold on;
ylabel('f [Hz]');

subplot(4,1,4);
plot(ts, delta * (180/pi));
hold on;
ylabel('delta [deg]');

xlabel('Time [s]');
