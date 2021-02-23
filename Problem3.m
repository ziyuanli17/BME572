%Problem 3
clear all;
load('spikes.mat')
load('stimulus.mat')

%Find S
S_bin = get_bin(0.0, 0.1, 0.1, 21);
sti_bool=[];
for b=S_bin.'
    sti_check = 0;
    for st = stimulus.'
        if b(2) >= st(1) && b(1) <= st(2)
            sti_check = 1;
            break;
        end
    end

    if sti_check == 1
        sti_bool=[sti_bool 1];
    else
        sti_bool=[sti_bool 0];
    end
end

S=[];
for t = get_bin(0, 20, 1, 201).'
    S=[S; sti_bool(t(1)+1:t(2)+1)];
end


% Get R
R=[];
[m, n] = size(spikes);

for t = get_bin(1.9, 2.0, 0.1, 20.1).'
    r_count=zeros(1,m);
    s_idx=1;   
    for trial=spikes.'
        for s=trial.'
            if s == 0
                break;
            end
            
            if t(2) >= s && t(1) <= s
                r_count(s_idx) = r_count(s_idx)+1;
            end
        end
        s_idx=s_idx+1;
    end
    R=[R; r_count];
end

%Find W
W1 = (pinv(S))*R(:, 1);
W2 = (pinv(S))*R(:, 2);
W3 = (pinv(S))*R(:, 3);
W4 = (pinv(S))*R(:, 4);
W=(W1+W2+W3+W4)./4;
prediction = NL_filter(S*W);

a=figure(1);
plot(1:21, flip(W));
xlabel('Time (s)');
ylabel('The Stimulus Weighting');
title('Linear Filter'); 
saveas(a, 'p3-1.png');

b=figure(2);
fplot(@(x) NL_filter(x));
xlabel('Input');
ylabel('Filter Resonse');
title('Non-linear Filter'); 
saveas(b, 'p3-2.png');

c=figure(3);
plot(1:182, R(:, 5));
xlabel('Time (s)');
ylabel('Spike Count');
title('Linear-Non-linear Spike Response'); 
hold on;
plot(1:182, prediction);
legend("Actual Response", "Predicted Response")
hold off;
saveas(c, 'p3-3.png');

% Generate the bin
function bin_arr = get_bin(init_start, init_end, incre, stop)
    bin_arr=[];
    while 1
        if init_end >= (stop + incre)
            break;
        end
        bin_arr = [bin_arr; [init_start init_end]];
        init_start=init_start+incre;
        init_end=init_end+incre;
    end
end


function arr_new=NL_filter(arr)
    for i=1:length(arr)
        if arr(i) < 0
            arr(i) = (-arr(i));
        end
        
        if arr(i) < 1
            arr(i) = arr(i)^3;
        end
    end
    arr_new=arr;
end


