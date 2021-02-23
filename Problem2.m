%Problem 2
clear all;
% Run the models aksed to implement
integrate_and_fire_neuron();
sinusoidal_IAF();
new_model();
sinusoidal_new_model();
two_neuron_oscillator();

%(a) Implement a simple integrate and fire neuron (IAF)
function integrate_and_fire_neuron
    %Initialize constants
    R=10;
    C=1;
    Vthr=5;
    Vspk=70;
    
    step_num = 100; % # of interations
    dt=1; %step size
    T = 0:step_num; %Time array

    % Initiliza array
    V = zeros(1,step_num);
    V(1)=0;
    
    I_arr=[0];
    skip=0;
    for t=1:(length(T)-1)
        %Get current
        if t >= 10 && t < 60
            I = 1;
        else
            I = 0;   
        end
        I_arr=[I_arr I];
        %Get V
        if skip == 1 %Skip if v was calculated (0)
            skip=0;
            continue;
        end

        if V(t)<Vthr %If V below threshold
            V(t+1) = V(t)+((I-(V(t)/R))/C)*dt;
        else %If V above threshold
            V(t+1) = Vspk;
            V(t+2) = 0;
            skip=1;
        end
    end
    
    %Plot the model
    a=figure(1);
    
    title('Integrate and Fire Neuron');
    yyaxis left;
    xlabel('Time (ms)');
    plot(T, V)
    ylabel('Voltage (mV)');
    
    yyaxis right;
    ylabel('Current (nA)');
    hold on
    plot(T, I_arr);
    legend('Voltage Response', 'Injected Current')
    saveas(a,'p2-1.png');
    hold off;
end

function sinusoidal_IAF
    %Initialize constants
    R=10;
    C=1;
    Vthr=5;
    Vspk=70;
    
    freq_arr = [1, 2, 5, 10, 20, 50, 100];
    spike_count_arr = zeros(1, length(freq_arr)); % For conting spikes at each freq
    
    step_num = 1000; % # of interations
    dt=1; %Time step size
    T = 0:step_num; %Time array

    % Initiliza array
    V = zeros(1,step_num);
    V(1)=0;
    
    skip=0;
    for f_idx = 1:(length(freq_arr))
        f = freq_arr(f_idx);
        I_arr=[0];
        for t=1:(length(T)-1)
            %Get sin current
            if t >= 0 && t < 1000
                I = sin(2*pi*f*(t*dt/1000));
            else
                I = 0;   
            end
            
            I_arr=[I_arr I];
            %Get V
            if skip == 1 %Skip if v was calculated (0)
                skip=0;
                continue;
            end

            if V(t)<Vthr %If V below threshold
                V(t+1) = V(t)+((I-(V(t)/R))/C)*dt;
                
            else %If V above threshold, spike
                V(t+1) = Vspk;
                V(t+2) = 0;
                % increase spike at freq by 1
                spike_count_arr(f_idx)=spike_count_arr(f_idx)+1;
                skip=1;
            end
        end
        
        %Plot V vs. t at a freq
        b=figure(f_idx+1); 
        
        title('IAF ' + string(freq_arr(f_idx)) + 'Hz');  
        yyaxis left;
        xlabel('Time (ms)');
        plot(T, V)
        ylabel('Voltage (mV)');
        hold on
        yyaxis right;
        plot(T, I_arr);
        ylabel('Current (nA)');
        legend('Voltage Response', 'Injected Current')
        saveas(b,string(f_idx)+'p2-2.png');
        hold off;
    end
    % Plot spike count vs. stimulus frequency
    c=figure(length(freq_arr)+2);
    plot(freq_arr, spike_count_arr)
    xlabel('Frequency (Hz)');
    ylabel('Spike Count');
    title('Spike Count vs. Stimulus Frequency');  
    saveas(c,'p2-3.png');
    
end

%Additional problems only for BME 572
%(b)
function new_model    
    step_num = 100; % # of interations
    dt=1; %Time step size
    T = 0:step_num; %Time array

    %Define constants
    a=0.02;
    b=0.2;
    c=-65;
    d=8;
    
    % Initiliza arrays
    v = zeros(1,step_num);
    u = zeros(1,step_num);
    v(1) =-65;
    u(1)=b*v(1);
    
    I_arr=[0];
    reset=0;
    for t=1:(length(T)-1)
        %Get sin current
        if t >= 10 && t < 60
            I = 10;
        else
            I = 0;   
        end
        I_arr=[I_arr I];
        %Reset
        if reset==1
            v(t+1) = c;
            u(t+1) = u(t)+d;
            % increase spike at freq by 1
            reset=0;
        else
            u(t+1) = u(t)+(a*(b*v(t)-u(t)))*dt;
            if v(t)>=30
                v(t)=30;
                v(t+1)=30;
                reset=1;
            else
                v(t+1) = v(t)+(0.04*(v(t)^2)+5*v(t)+140-u(t)+I)*dt;
            end
        end
    end 
    
    %Plot V vs. t
    d=figure(9);
    xlabel('Time (ms)');
    yyaxis left;
    plot(T, v)
    ylabel('Voltage (mV)')
    hold on
    plot(T, u);
    yyaxis right;
    ylabel('Current (nA)');
    plot(T, I_arr);
    
    legend("v", "u", 'Injected Current')
    title('New Model Neuron');  
    saveas(d,'p2-4.png');
    hold off;   
end

function sinusoidal_new_model    
    step_num = 1000; % # of interations
    dt=1; %Time step size
    T = 0:step_num-1; %Time array

    %Define constants
    a=0.02;
    b=0.2;
    c=-65;
    d=8;
    
    % Initilize arrays
    v = zeros(1,step_num);
    u = zeros(1,step_num);
    v(1) =-65;
    u(1)=b*v(1);
    
    freq_arr = [1, 2, 5, 10, 20, 50, 100];
    spike_count_arr = zeros(1, length(freq_arr)); % For conting spikes at each freq
    
    reset=0;
    for f_idx = 1:(length(freq_arr))
        f = freq_arr(f_idx);
        I_arr=[0];
        for t=1:(length(T)-1)
            %Get sin current
            if t >= 0 && t < 1000
                I = 10*sin(2*pi*f*(t*dt/1000));
            else
                I = 0;   
            end
            I_arr=[I_arr I];
            
            %Reset
            if reset==1
                v(t+1) = c;
                u(t+1) = u(t)+d;
                % increase spike at freq by 1
           
                reset=0;
                
            else
                u(t+1) = u(t)+(a*(b*v(t)-u(t)))*dt;
                if v(t)>=30
                    spike_count_arr(f_idx)=spike_count_arr(f_idx)+1;
                    
                    v(t)=30;
                    v(t+1)=30;
                    reset=1;
                else
                    v(t+1) = v(t)+(0.04*(v(t)^2)+5*v(t)+140-u(t)+I)*dt;
                end
            end
        end
        
        %Plot V vs. t at a freq
        e=figure(length(freq_arr)+f_idx+3);
        xlabel('Time (ms)');
        yyaxis left;
        plot(T, v)
        ylabel('Voltage (mV)');
        hold on
        plot(T, u);
        yyaxis right;
        plot(T, I_arr);
        ylabel('Current (nA)');
        legend("v", "u", 'Injected Current')
        title('New Model ' + string(freq_arr(f_idx)) + 'Hz');   
        hold off;
        saveas(e,string(f_idx)+'p2-5.png');
    end
    
    % Plot spike count vs. stimulus frequency
    f=figure(2*length(freq_arr)+4);
    plot(freq_arr, spike_count_arr)
    xlabel('Frequency (Hz)');
    ylabel('Spike Count');
    title('Spike Count vs. Stimulus Frequency'); 
    saveas(f,string(f_idx)+'p2-6.png');
end


%(c)Construct a two-neuron oscillator using reciprocal inhibition
function two_neuron_oscillator
    step_num = 1500; % # of interations
    dt=1; %Time step size
    T = 0:step_num; %Time array
    
    %Define constants
    C=1;
    R=10;
    Vrest=0;
    Vspk=70;
    tauthre=50;
    Einh=-15;
    tausyn=15;
    gpeak=0.1;
    
    %Initialize an array of structs for 2 neurons
    neuron_arr=[[] []];
    for i=1:2
        neuron_arr(i).v = zeros(1,step_num);
        neuron_arr(i).v(1) = Einh;
        neuron_arr(i).theta = zeros(1,step_num);
        neuron_arr(i).theta(1) = Vrest;
        neuron_arr(i).z = zeros(1,step_num);
        neuron_arr(i).z(1) = 0.2;
        neuron_arr(i).g = zeros(1,step_num);
        neuron_arr(i).g(1) = 0.2;
    end
    
    reset_arr = [0 0];
    for t=1:(length(T)-1)
        %Get currents
        if t >= 0 && t < 1500
            I1=1.1;
            I2=0.9;
        else
            I1=0;
            I2=0;  
        end
        
        for i=1:(length(neuron_arr))
            j=2;
            I=I1;
            if i == 2
                I=I2;
                j = 1;
            end
            
            % Euler's method
            neuron_arr(i).v(t+1) = neuron_arr(i).v(t) + (((-neuron_arr(i).v(t)/R) - neuron_arr(i).g(t)*(neuron_arr(i).v(t) - Einh) + I)/C)*dt;                        
            neuron_arr(i).theta(t+1) = neuron_arr(i).theta(t) + ((-neuron_arr(i).theta(t) + neuron_arr(i).v(t))/tauthre)*dt;
            neuron_arr(i).z(t+1) = neuron_arr(i).z(t) + (-neuron_arr(i).z(t)/tausyn) + (gpeak/(tausyn/exp(1)))*(neuron_arr(j).v(t)==Vspk)*dt;
            neuron_arr(i).g(t+1) = neuron_arr(i).g(t) + ((-neuron_arr(i).g(t)/tausyn) + neuron_arr(i).z(t))*dt;      
            
            if reset_arr(i) == 1
                neuron_arr(i).v(t+1) = Einh;
                reset_arr(i)=0;
            end
            
            if neuron_arr(i).v(t+1) >= neuron_arr(i).theta(t+1) %when fires
                neuron_arr(i).v(t+1) = Vspk;
                reset_arr(i) = 1;
            end
        end
    end 
    
    %Plot V vs. t
    g=figure(19);
    plot(T, neuron_arr(1).v)       
    xlabel('Time (ms)');
    ylabel('Voltage (mV)');
    title('Two Neuron Oscillator');  
    hold on;
    plot(T, neuron_arr(2).v);
    legend("Neuron 1", "Neuron 2")
    hold off; 
    saveas(g, 'p2-7.png');
end

