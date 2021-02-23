%Problem 4
lambda_max = max(prediction);
spike_trains = zeros(length(prediction),100);
for t = 1:1000
    spike_trains(:,t) = rand(length(prediction),1) <= (prediction)/lambda_max;
end

seq=[0:0.1:18.1];
% Plot 4 spike trains
a=figure(1);
plot(seq, spike_trains(:,1));
xlabel('Time (s)');
ylabel('Spike Occurence');
title('Spike Train 1'); 
saveas(a, 'p4-1.png');

b=figure(2);
plot(seq, spike_trains(:,2));
xlabel('Time (s)');
ylabel('Spike Occurence');
title('Spike Train 2'); 
saveas(b, 'p4-2.png');

c=figure(3);
plot(seq, spike_trains(:,3));
xlabel('Time (s)');
ylabel('Spike Occurence');
title('Spike Train 3'); 
saveas(c, 'p4-3.png');

d=figure(4);
plot(seq, spike_trains(:,3));
xlabel('Time (s)');
ylabel('Spike Occurence 4');
title('Spike Train 4'); 
saveas(d, 'p4-4.png');

% Get ISI
ISI=[];
for t = 1:1000
    spike_col = spike_trains(:,t)';
    t_diff=[];
    interval_count = 0;
    t_pre = 0;
    for s = spike_col
       interval_count=interval_count+1;
       if s == 1
           t_diff=[t_diff (interval_count-t_pre)];
           t_pre=interval_count;
       end
    end   
    ISI=[ISI t_diff];
end

%ISI distribution
e=figure(5);
histogram(ISI, 50);
xlabel('Interspike Interval (S)');
ylabel('Interval Difference');
title('ISI distribution'); 
saveas(e, 'p4-5.png');

%Fano
fano_factor = var(spike_trains)/mean(spike_trains);
disp("The Fano Factor is: " + string(fano_factor));

%Coefficient of Variation
coefficient_of_variation = std(ISI)/mean(ISI);
disp("Coefficient of Variation for ISI is: " + string(coefficient_of_variation));




