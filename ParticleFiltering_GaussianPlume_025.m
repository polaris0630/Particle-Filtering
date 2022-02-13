
% Particle filtering

p = 2; % initial state
Q = 0.001; % process noise covariance
R = 0.001; % measurement noise covariance
SL = 500; % simulation length
 
N = 200; % number of particles in the particle filter

x = 2500; % distance between source and city
h = 100; % height of dispersion source
miu = 0.25; % average wind speed
s = 0.3; % original dispersion speed

phat = p;
phatPart = p;

lamda = 1.527 * 10^(-6); % decay constant
sigmap = 0.01; % standard deviation of concentration of source

sigmax = (0.32 * x) / (1 + 0.0004 * x)^0.5;
sigmaz = 0.24 * x * (1 + 0.001 * x)^0.5;

t = 0; % time

% Initialize the particle filter.
for i = 1 : N
    ppart(i) = p + sqrt(sigmap) * randn; % randn: generate random numbers which follow N(1,2^0.5).
end

pArr = [];
CArr = [0];
phatArr = [p];
phatPartArr = [];

pArract = [];

close all;
 
floor((x/(s+miu))/3600)

time = floor((x / (s + miu)) / 3600);

for t = 1 : SL + time
    
    
    if t <= time
            
    % Non-linear State Space Model Function
    p = p;
    C = 0;
   
    pArr = [pArr p];
    CArr = [CArr C];
    

    else
        
    factor = (sigmax^2 * h^2) / (sigmaz^2 * x^2 + sigmax^2 * h^2);
    
    p = p / exp(lamda * 3600) + sqrt(Q) * randn  ;
    C = pArr(t - time) * factor * exp(-(x + miu * 3600 * t)^2 / (2 * sigmax^2)) + sqrt(R) * randn;
   
    % Particle filter
    for i = 1 : N
        ppartminus(i) = ppart(i) / exp(lamda * 3600) + sqrt(Q) * randn;
        Cpart = ppartminus(i) * factor * exp(-(x + miu * 3600 * t)^2 / (2 * sigmax^2)) + sqrt(R) * randn;
        vhat = Cpart - C;
        q(i) = (1 / (sqrt(R) * sqrt(2*pi))) * exp(-vhat^2 / (2 * R));
    end
    
    % Normalize the weight for each estimate.
    qsum = sum(q);
    
    for i = 1 : N
        q(i) = q(i) / qsum;
    end
    
    % Resampling
    for i = 1 : N
        u = rand; % uniform random number between 0 and 1
        qtempsum = 0;
        for j = 1 : N
            qtempsum = qtempsum + q(j);
            if qtempsum >= u
                ppart(i) = ppartminus(j);
                qArr(i) = q(j);
                break;
            end
        end
    end
    
    % The particle filtering estimation is the weighted mean of the particles.
    phatPart = sum(ppart.*qArr) / sum(qArr);
    
    % Save data in arrays for later plotting
    pArr = [pArr p];
    CArr = [CArr C];
    phatArr = [phatArr phat];
    phatPartArr = [phatPartArr phatPart];

    
    end;
  
    
end

for l = 1 : SL 
    
pArract = [pArract pArr(l)];

end

k = 1 : SL;

figure;
plot(k, pArract, 'b.', k, phatPartArr, 'k-');
set(gca,'FontSize',12); set(gcf,'Color','White'); 
xlabel('Time Step'); ylabel('State');
legend('True State', 'Particle Filtering Estimation'); 

pArract = transpose(pArract);
phatPartArr = transpose(phatPartArr);

xlswrite('0.25real.xls', pArract);
xlswrite('0.25estimate.xls', phatPartArr);