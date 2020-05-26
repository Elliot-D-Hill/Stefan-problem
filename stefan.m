close all; clear; clc;

% Stefan moving boundary value problem

% choose interpolation strategy from: 'taylor' or 'poly'
interpStep = 'poly';
% choose extrapolation strategy from: 'taylor', 'taylorSafe' or 'poly'
extrapStep = 'taylorSafe';

% Initialize Data
a = 0; % start of x domain
beta = 0.01; % velocity of h(t)
TM = 0;
L = 2; % length of rod
tEnd = 500; % end of time interval
nu = 0.01; % heat constant
N = 500; % number of subintervals in space
dx = L/N;
h = 4.123*dx; % initial depth of ice
alpha = 0.5;
dt = alpha * dx^2 / nu;
M = fix(tEnd / dt); % number of time steps
mGrid = fix(h / dx);

x = a:dx:L; % grid of x_i

% exact solution parameters
q1 = 0.0901894;
t0 = ((h/2/q1)^2) / nu;
htExact = NaN(M,1);

T = zeros(N+1, 1); % intial conditions
f = @(x) exp(-x) - exp(-h); % intial conditions

T(1:mGrid) = f(x(1:mGrid));
TA = f(0);
T(1) = TA;

t = NaN(M,1);
ht = NaN(M,1);

% updating T with heat equation
for m = 1:M % for each time step
    T(1)=TA;
    for n = 2:mGrid % for each x in grid
        T(n) = T(n) + alpha * (T(n+1) - 2 * T(n) + T(n-1));
    end
    
    % distance between h^n and mth grid point
    eta = h - (mGrid * dx); 
    
    switch extrapStep
        
        case 'taylor'
            Tx = gradExterp(eta, mGrid, dx ,TM, T);
        case 'taylorSafe'
            % guard against small eta
            if eta < 0.1
                Tx = (T(mGrid) - T(mGrid-1)) / dx;
            else
                Tx = gradExterp(eta, mGrid, dx ,TM, T);
            end 
        case 'poly'
            % find gradient by quadratic polynomial extrapolation
            Tx = ((1 + 2*eta)*T(mGrid-2) - 4*T(mGrid-1)*(1+eta) + T(mGrid)*(3+2*eta)) / (2*dx);
    end
    
    % update h^n
    h = h - (dt * beta * Tx);
    
    % check to see if h^n crossed a grid point
    if h > (mGrid + 1) * dx
        % interpolate to find temp at new grid point
        switch interpStep
            case 'taylor'
                % Taylor series approximation
                eta = h - ((mGrid + 1) * dx);
                % coefficients for interpolation
                a1 = 2 / (eta + 2) * (eta + 1);
                a2 = (2 * eta) / (eta + 1);
                a3 = -eta / (eta + 2);
                
                T(mGrid+1) = a1 * TM + a2 * T(mGrid) + a3 * T(mGrid-1);
            case 'poly'
                % polynomial interpolation
                T(mGrid+1) = 3*T(mGrid) - 3*T(mGrid-1) + T(mGrid-2);
        end
        
        mGrid = mGrid + 1;
    end
    
    if h > L - dx % if close to end of grid
        break
    end
    
    t(m) = m*dt;
    ht(m) = h;
    htExact(m) = 2*q1*sqrt(nu*(t(m) + t0));
end

figure
% exact solution should look like square root
h1 = plot(t, ht, "LineWidth", 2);
set(h1, 'markerfacecolor', get(h1, 'color'));
hold on
h2 = plot(t, htExact, "LineWidth", 2, 'LineStyle','--');
set(h2, 'markerfacecolor', get(h2, 'color'));
legend("Numerical", "Exact")

xlabel('t', 'FontSize', 16)
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')
ylabel('h(t)', 'FontSize', 16)

% create a new pair of axes inside current figure
axes('position',[.5 .175 .35 .35])
box on % put box around new pair of axes
indexOfInterest = (t < 11*pi/8) & (t > 9*pi/8); % range of t near perturbation
plot(t(indexOfInterest),ht(indexOfInterest), "LineWidth", 2) % plot on new axes
axis tight

function Tx = gradExterp(eta, mGrid, dx ,TM, T)
% coefficients for extrapolation
c1 = ((2 * eta) + 1) / (eta * dx * (eta+1));
c2 = -(eta + 1) / (eta * dx);
c3 = eta / ((eta * dx) + dx);

% extrapolation step to find unknown gradient
Tx = c1 * TM + c2 * T(mGrid) + c3 * T(mGrid-1);
end


