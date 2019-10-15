%%%%%%%%%%%%%%%%%%%%%%%% DECLARATION OF VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%

m = 100;            % total number of nodes
lambda = 1/exp(1);  % overall arrival rate in [packets/slot]
T = 1000;           % simulation duration in [slots]

q_r = 0.01;                 % probability that a backlogged node retransmits in the next slot (user-defined param)
q_a = 1 - exp(-lambda/m);   % probability that an unbacklogged node transmits in the next slot

% Auxiliar counters:
tx_attemps_counter = 0;
successes_counter = 0;
collisions_counter = 0;

nodes_state = zeros(1,m);
% 0 = idle node waiting to receive new packets
% 1 = unbacklogged node ready to tx new packet
% 2 = backlogged node ready to retx packet
% 3 = backlogged node waiting

num_of_backlogged = zeros(1,T); % to plot backlogged VS slot and histogram
num_of_arrivals = zeros(1,T); % to plot packets_in VS slot
num_of_departures = zeros(1,T); % to plot packets_out VS slot

%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

prompt = 'Do you use the default values for arrival rate and retx probability? Y/N [Y]: ';
str = input(prompt,'s');
if isempty(str)
    str = 'Y';
end

if str~='Y' && str~='y'
    lambda = input('Arrival rate: ');
    q_r = input('Retx probability: ');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s = 1; % current slot
while s < T
    % Iterate over all nodes to detect new arrivals:
    for node = 1:m
        if nodes_state(1,node)==0 && rand(1)<=q_a
            % The idle node under study has received a new packet to tx
            nodes_state(1,node) = 1; % update node's state to "ready_to_tx"
            num_of_arrivals(1,s) = num_of_arrivals(1,s) + 1; % count new arrival
        elseif nodes_state(1,node)==3 && rand(1)<=q_r
            % The backlogged node under study is ready to retx (waiting time concluded)
            nodes_state(1,node) = 2; % update node's state to "ready_to_retx"
        end
    end
    
    % All nodes with state 1 or 2 will transmit in the current slot:
    tx_attempts_current_slot = sum(nodes_state==1) + sum(nodes_state==2);
    tx_attemps_counter = tx_attemps_counter + tx_attempts_current_slot;
    
    successful_retx = 0; 
    [~,tx_ids] = find(nodes_state==1);
    [~,retx_ids] = find(nodes_state==2);
    if tx_attempts_current_slot == 1 % if only one packet transmits...
        successes_counter = successes_counter + 1;
        num_of_departures(1,s) = num_of_departures(1,s) + 1;
        % The only node that transmitted becomes idle, update state:
        if isempty(retx_ids) % successful tx
            nodes_state(1,tx_ids) = nodes_state(1,tx_ids) - 1;
            num_of_backlogged(1,s+1) = num_of_backlogged(1,s);
        else % successful retx
            nodes_state(1,retx_ids) = nodes_state(1,retx_ids) - 2;
            num_of_backlogged(1,s+1) = num_of_backlogged(1,s) - 1;
            successful_retx = 1;
        end
    elseif tx_attempts_current_slot > 1 % if several packets transmit...
        collisions_counter = collisions_counter + 1;
        % Only those nodes that were unbacklogged become backlogged:
        num_of_backlogged(1,s+1) = num_of_backlogged(1,s) + length(tx_ids);
        % All nodes involved in the collision change to "waiting" state:
        nodes_state(1,tx_ids) = 3;
        nodes_state(1,retx_ids) = 3;
    else % if no transmission occurs...
        num_of_backlogged(1,s+1) = num_of_backlogged(1,s);
    end
    
    num_of_arrivals(1,s+1) = num_of_arrivals(1,s);
    num_of_departures(1,s+1) = num_of_departures(1,s);
    
    s = s + 1;
end

% Simulation results:
traffic_offered = tx_attemps_counter / T;
collision_prob = collisions_counter / T;
freq_of_success = successes_counter / T; % throughput
fprintf('Collision Probability: %d\nTraffic Offered: %d\n',collision_prob,traffic_offered);

% Calculate the steady-state probabilities of the Markov chain, as the backlog frequencies
[counts,~] = histcounts(num_of_backlogged);
steady_state_probs = counts / sum(counts);

% Use the theoretical formulas to calculate G and Ps as a function of the backlog:
aux = 1:m;
G = (m*ones(1,m)-aux)*q_a+aux*q_r; % attempt rate
Ps = G.*exp(-G); % probability of success

% Combine all the results to calculate the average probability of success.
steady_state_probs = [steady_state_probs, zeros(1,m-length(steady_state_probs))];
avg_prob_success = dot(steady_state_probs,Ps);

% Compare it with the value derived directly from the simulation results, as frequency of success.
fprintf('\nFrequency of success: %d\n',freq_of_success);
fprintf('Average probability of success: %d\n',avg_prob_success);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
h = histogram(num_of_backlogged);
xlabel('Num. backlogged nodes')
ylabel('Counts')
title('histogram of the backlog')

figure(2)
ex3g1 = subplot(2,1,1);
plot(ex3g1,1:T,num_of_backlogged)
xlabel(ex3g1,'Slot Number')
ylabel(ex3g1,'Num. backlogged nodes')
title(ex3g1,'system''s backlog VS slot number')
ex3g2 = subplot(2,1,2);
plot(ex3g2,1:T,num_of_arrivals)
hold on
plot(ex3g2,1:T,num_of_departures)
xlabel(ex3g2,'Slot Number')
ylabel(ex3g2,'Num. of arrivals/departures')
hold off
title(ex3g2,'Arrivals and Departures over time')
legend(ex3g2,'Number of Arrivals','Number of Departures')
legend('Location','northwest')