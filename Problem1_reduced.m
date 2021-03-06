%Problem 1. (a) 
%Implement the HH model in Matlab.

step_num = 100000;

% Initiliza arrays
n = zeros(1,step_num);
m = zeros(1,step_num);
h = zeros(1,step_num);
V = zeros(1,step_num);

ConducK = zeros(1,step_num+1);
ConducNa = zeros(1,step_num+1);
I_k = zeros(1,step_num+1);
I_Na = zeros(1,step_num+1);
I_leak = zeros(1,step_num+1);
I_cap = zeros(1,step_num+1);

% Define constants
Cm = 1;
I = 0;
Vl = -61;
Vk = -77;
VNa = 55;
gL = 0.3;
gK = 36;
gNa = 120;

dt=0.001; % step's size
T = 0:100000;

% Define init values
n(1) = 0.3;
m(1) = 0.06;
h(1) = 0.6;
V(1) = -65;
ConducK(1) = (gK*(n(1)^4));
ConducNa(1) = (gNa*(m(1)^3)*h(1));

I_k(1) = ConducK(1) * (V(1)-Vk);
I_Na(1) = ConducNa(1) * (V(1)-VNa);
I_leak(1) = gL * (V(1)-Vl);




for t=1:(length(T)-1)
    if t==60000
        I=30;
    end
    
    if t==63000
        I=0;
    end
    
    %Define alpha and betas
    an = (0.01 * ((-65-V(t))+10))/(exp(((-65-V(t))+10)/10)-1);
    bn = 0.125*exp((-65-V(t))/80);

    am = 0.1 * ((-65-V(t)) + 25)/(exp(((-65-V(t))+25)/10)-1);
    bm = 4*exp((-65-V(t))/18);
  
    ah = 0.07 * exp((-65-V(t))/20);
    bh = 1/(exp(((-65-V(t))+30/10))+1);
    
    % Conductances
    ConducK(t) = gK*(n(t)^4);
    ConducNa(t) = gNa*(m(t)^3)*(0.89-(1.1*n(t)));
    
    % Currents
    I_k(t) = ConducK(t) * (V(t)-Vk);
    I_Na(t) = ConducNa(t) * (V(t)-VNa);
    I_leak(t) = gL * (V(t)-Vl);
    if t>1
        I_cap(t) = (V(t)-V(t-1))/0.001;
    end
    
    % Membrane potential
    V(t+1) = V(t)+ ((I - I_k(t) - I_Na(t) - I_leak(t))/Cm)*dt;

    % HH variables
    n(t+1) = n(t) + (an*(1-n(t)) - bn*n(t))*dt;
    m(t+1) = am/(am+bm);
end

a=figure(1);
plot(T*0.001,V);
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
title('Action potential');
saveas(a,'p1-3-1.png');

b=figure(2);
plot(T*0.001,m);
xlabel('Time (ms)');
ylabel('HH Variables');
title('Hodgkin-Huxley Variables');  
hold on
plot(T*0.001,n);
legend("m", "n");
hold off
saveas(b,'p1-3-2.png');

%Current
c=figure(3);
disp(length(T))
disp(length(I_k))
plot(T*0.001,I_k);
title('Membrane Currents');
xlabel('Time (ms)');
ylabel('Current (mA/cm^2)');
hold on;
plot(T*0.001,I_Na);
plot(T*0.001,I_leak);
plot(T*0.001,I_cap);
legend("I_k", 'I_N_a', "I_l_e_a_k", "I_c_a_p");
hold off;
saveas(c,'p1-3-3.png');

%Conductance
d=figure(4);
plot(T*0.001,ConducK);
title('Ion Comductance');
xlabel('Time (ms)');
ylabel('Conductance (mS/cm^2)');
hold on;
plot(T*0.001,ConducNa);
legend("G_K", "G_N_a")
hold off;
saveas(d,'p1-3-4.png');