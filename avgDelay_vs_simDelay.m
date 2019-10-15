%%%%%%%%%%%%%%%%%%%%%%%% DECLARATION OF VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%

m = 100;            % total number of nodes
q_r = 1;            % probability that a backlogged node retransmits in the next slot (user-defined param)
T = 1000;           % simulation duration in [slots]
feedback = 0;       % 0=idle , 1=successful_tx , 2=collision

tx_attemps_counter = 0;
successes_counter = 0;
collisions_counter = 0;

nodes_state = zeros(1,m);
n_hats = ones(1,T);
last_tx_slot = zeros(1,m);
success_packets_delay = zeros(1,T);
num_of_backlogged = zeros(1,T); % to plot backlogged VS slot and histogram
num_of_arrivals = zeros(1,T); % to plot packets_in VS slot
num_of_departures = zeros(1,T); % to plot packets_out VS slot

arr_rate = 0.05:0.05:0.35;
approx_avg_delay = zeros(1,length(arr_rate));
real_avg_delay = zeros(1,length(arr_rate));
i = 1;
while i<=length(arr_rate)
    approx_avg_delay(1,i) = (exp(1)-0.5)/(1-arr_rate(i)*exp(1)) - ((exp(1)-1)*(exp(arr_rate(i))-1))/(arr_rate(i)*(1-(exp(1)-1)*(exp(arr_rate(i))-1)));
    q_a = 1 - exp(-arr_rate(i)/m);
    
    s = 1; % current slot
    while s < T
        % Iterate over all nodes to detect new arrivals:
        for node = 1:m
            if nodes_state(1,node)==0 && rand(1)<=q_a
                % The idle node under study has received a new packet to tx
                num_of_arrivals(1,s) = num_of_arrivals(1,s) + 1; % count new arrival
                if rand(1)<=q_r
                    % The new arrival is transmitted in the next slot
                    nodes_state(1,node) = 1; % update node's state to "ready_to_tx"
                    last_tx_slot(1,node) = s; % keep time of first attempt to calculate delay later
                else
                    % The new arrival will have to wait to be transmitted
                    nodes_state(1,node) = 4; % update node's state to "waiting"
                end
            elseif nodes_state(1,node)==3 && rand(1)<=q_r
                % The backlogged node under study is ready to retx (waiting time concluded)
                nodes_state(1,node) = 2; % update node's state to "ready_to_retx"
            elseif nodes_state(1,node)==4 && rand(1)<=q_r
                % The backlogged node under study is ready to retx (waiting time concluded)
                nodes_state(1,node) = 1; % update node's state to "ready_to_tx"
            end
        end
        
        % All nodes with state 1 or 2 will transmit in the current slot:
        tx_attempts_current_slot = sum(nodes_state==1) + sum(nodes_state==2);
        tx_attemps_counter = tx_attemps_counter + tx_attempts_current_slot;
        
        [~,tx_ids] = find(nodes_state==1);
        [~,retx_ids] = find(nodes_state==2);
        if tx_attempts_current_slot == 1 % if only one packet transmits...
            successes_counter = successes_counter + 1;
            feedback = 1; % indicates successful departure
            num_of_departures(1,s) = num_of_departures(1,s) + 1;
            % The only node that transmitted becomes idle, update state:
            if isempty(retx_ids) % successful tx
                nodes_state(1,tx_ids) = 0;
                success_packets_delay(successes_counter) = s - last_tx_slot(tx_ids);
                num_of_backlogged(1,s+1) = num_of_backlogged(1,s);
            else % successful retx
                nodes_state(1,retx_ids) = 0;
                success_packets_delay(successes_counter) = s - last_tx_slot(retx_ids);
                num_of_backlogged(1,s+1) = num_of_backlogged(1,s) - 1;
            end
        elseif tx_attempts_current_slot > 1 % if several packets transmit...
            collisions_counter = collisions_counter + 1;
            feedback = 2; % indicates collision
            % Only those nodes that were unbacklogged become backlogged:
            num_of_backlogged(1,s+1) = num_of_backlogged(1,s) + length(tx_ids);
            % All nodes involved in the collision change to "waiting" state:
            nodes_state(1,tx_ids) = 3;
            nodes_state(1,retx_ids) = 3;
        else % if no transmission occurs...
            feedback = 0; % indicates idle slot
            num_of_backlogged(1,s+1) = num_of_backlogged(1,s);
        end
        
        num_of_arrivals(1,s+1) = num_of_arrivals(1,s);
        num_of_departures(1,s+1) = num_of_departures(1,s);
        
        s = s + 1;
        if feedback==2
            n_hats(1,s) = n_hats(1,s-1) + arr_rate(i) + 1/(exp(1)-2);
        else
            n_hats(1,s) = max(arr_rate(i),n_hats(1,s-1)+arr_rate(i)-1);
        end
        q_r = min(1,1/n_hats(1,s));
    end
    last_tx_slot = zeros(1,m);
    real_avg_delay(1,i) = mean(success_packets_delay);
    i = i + 1;
end

% Plot the average delay for a successful packet tx VS the arrival rate
% using the formula obtained by the approximate delay analysis of the Pseudo-Bayesian algorithm.
% Plot as well the average delay for a successful packet tx,
% calculated from the simulation results of stabilized ALOHA.
figure(1)
plot(arr_rate,approx_avg_delay)
hold on
plot(arr_rate,real_avg_delay)
xlabel('Arrival Rate (lambda)')
ylabel('Average Delay')
hold off
title('Average delay for a successful packet tx')
legend('From theoretical formula','From simulation results')
legend('Location','northwest')