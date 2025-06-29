clc
clear all
close all

settings = calibrationSettings %load settings from different .m file
empVolatilitydata = load("empVolatilitySurfaceData.mat"); %load the providid volatility matrix
empVolatilitydata = empVolatilitydata.data %take only the data part of the mat file

%% params

tradingDaysInYear = settings.tradingDaysInYear;
n = settings.n; %setting for the FFT pricing algorithm
model = settings.model; %Name of the pricing al

%Initialize the model parameters, these could be "yesterdays" values or
%random values within the range of possible values
kappa = settings.parameters0(1);
theta = settings.parameters0(2);
eta = settings.parameters0(3);
rho = settings.parameters0(4);
V0 = settings.parameters0(5);

%put the initialized values in an array for later
param0 = [kappa theta eta rho V0];

%Initialize the minimiums and maximiums of paramter values
minKappa = settings.minKappa;
maxKappa = settings.maxKappa;
minTheta = settings.minTheta;
maxTheta = settings.maxTheta;
minEta = settings.minEta;
maxEta = settings.maxEta;
minRho = settings.minRho;
maxRho = settings.maxRho;
minV0 = settings.minV0;
maxV0 = settings.maxV0;


OptSettings = settings.calibrOptions; %settings for the optimization algorithm

r = empVolatilitydata.r; %Risk free interst rate
S0 = empVolatilitydata.S0; %price of the underlying asset now
K = empVolatilitydata.K; %strike prices
T = empVolatilitydata.T; % Times to maturities

IVolSurf = empVolatilitydata.IVolSurf; %volatilities at different K/T pairs



%% Optimization
%set model parameter limits into arrays
var_max_limits = [maxKappa, maxTheta,maxEta, maxRho, maxV0]; 
var_min_limits = [minKappa, minTheta,minEta, minRho, minV0];

%Control the output of text and graphs during optimization
settings.standardErrors = false;
settings.indPlotSurface = 1;

%define the loss function that optimizes the mean square error between the
%"observed volatility surface and the model output volatility surface
fun = @(x) lossFunction1(x, var_min_limits, var_max_limits,  model, n, S0, K, T, r, IVolSurf, settings);
[param_final,fval,exitflag,output] = patternsearch(fun, param0)
%[param_final, fFinal, exitFlag, output, grad, hessian] = fminunc(fun, param0, OptSettings); %Matlab function to optimize the loss value of a function using gradients

% Results
disp(['Starting values: ', num2str(param0)]);
disp(['Optimized values: ', num2str(param_final)]);
%% Statistical evaluation of parameters (standard error)
%I dont want to overwrite n from previous the amount of elements in the
%volatility surface matrix is called n2
n2 = numel(IVolSurf);
f = fun(param_final)*n2; %make the mean standard error into standard error

%Define the loss function for the calculation of standard errors for the parameters. 
% In practice this gives out the model volatility surface instead of a loss value. 
settings.standardErrors = true;
settings.indPlotSurface = 0; %show the surfaces when calling the funtion
fis = @(x) lossFunction1(x, var_min_limits, var_max_limits,  model, n, S0, K, T, r, IVolSurf, settings);

p = length(param0); %how many params there are
J = jacobian(fis, param_final); %Calculate jacobian matrix with jacobian.m file

%calculating standard errors with jacobian matrix
sigma2 = f/(n2-p);
Sigma = sigma2*inv(J'*J);
se = sqrt(diag(Sigma))';

%Displaing the error values
disp(['Standard error values: ', num2str(se)]);




%% Price an exotic option with the optimized parameter values

%This is a down-and-in arithmetic Asian average strike call option

M = 100000; %number of simulations
H = S0*0.85; %Barrier that has to be hit in order to get payoff
T2 = 1-0; %Time to maturity (T-t)
dt = 1/tradingDaysInYear; %time intervals to calculate prices for (1 day)
N = T2/dt; %Number time points calculated within 1 simulation. Now its 252 same as tradingDaysInYear

%Get the optimized parameters out of the array for clarity
kappa_final = param_final(1);
theta_final = param_final(2);
eta_final = param_final(3);
rho_final = param_final(4);
V0_final = param_final(5);

gamma = 0; % Assume risk-neutrality so gamma is 0

%initialize arrays for payoffs, antithetic payoffs and the array for
%combination of them
payoffs = [];
payoffs_a = [];
W = [];

 
for sim = 1:M
    %Initialize stock price arrays for each simulation and their starting
    %prices
    S = [];
    S(1) = S0;
    %same for antithetic prices
    S_a = [];
    S_a(1) = S0;
    
    % Initialize the squared instantaneous volatility
    V = [];
    V(1) = V0_final;
    %same for the antithetic version
    V_a = [];
    V_a(1) = V0_final;
    
    % Simulate the path for stock prices and the squaared instantaneous
    % volatility, because of the barrier condition
    for t = 2:N 
        %create the random variable necessery for the Milstein scheme of
        %Heston
        e12 = randn(); %Initialize independent normally distributed random variable
        e12_a = -e12; %antithetic version of previous

        epsilon1 = randn(); %Initialize another independent normally distributed random variable
        epsilon1_a = -epsilon1; %antithetic version of previous

        %Correlated variables using correlation coefficent rho with the
        %epsilon1 and epsilon1_a
        epsilon2 = rho_final*epsilon1 + sqrt(1-rho_final^2)*e12;
        epsilon2_a = rho_final*epsilon1_a+sqrt(1-rho_final^2)*e12_a;

        %Calculate and store squared instantaneous volatility for this time
        %step
        V(t) = max(V(t-1)+kappa_final*(theta_final-V(t-1))*dt+eta_final*sqrt(V(t-1)*dt)*epsilon2+0.25*eta_final^2*dt*(epsilon2^2-1),0);
       
        %Calculate and store the stock price for this time step
        S(t) = S(t-1)*exp((r+gamma*V(t-1)-0.5*V(t-1))*dt+sqrt(V(t-1))*sqrt(dt)*epsilon1);

        %Calculate and store the same for the antithetic versions
        V_a(t) = max(V_a(t-1)+kappa_final*(theta_final-V_a(t-1))*dt+eta_final*sqrt(V_a(t-1)*dt)*epsilon2_a+0.25*eta_final^2*dt*(epsilon2_a^2-1),0);
        S_a(t) = S_a(t-1)*exp((r+gamma*V_a(t-1)-0.5*V_a(t-1))*dt+sqrt(V_a(t-1))*sqrt(dt)*epsilon1_a);

    end
    
    %Calculate average and minimum of the stock price
    S_avg = mean(S);
    S_min = min(S(:));

    %Check if barrier was hit during the simulation if so: calculate
    %payoff, if not payoff is 0
    if S_min < H
        payoffs(sim) = exp(-r*T2)*max(S(end)-S_avg,0); %calculate the arithmetic Asian average strike call option and
        % Convert the payoff to current prices
    else
        payoffs(sim) = 0;
    end

    %Calculate average and minimum of the antithetic stock price
    S_a_avg = mean(S_a);
    S_a_min = min(S_a);
    %Check if barrier was hit during the antithetic simulation if so: calculate
    %payoff, if not payoff is 0
    if S_a_min < H
        payoffs_a(sim) = exp(-r*T2)*max(S_a(end)-S_a_avg,0); %calculate the arithmetic Asian average strike call option and
        % Convert the payoff to current prices
    else
        payoffs_a(sim) = 0;
    end
    W(sim) = 0.5*(payoffs(sim)+payoffs_a(sim)); %calculate the antithetic option payoff
end

price = mean(W(:)) % calculat the average antithetic option payoff over all the simulations

 

 

 

 

%% Loss function

function loss = lossFunction1(param, var_min, var_max, model, n, S0, K, T, r, IVolSurf, settings)

    % Get the parameters out of the array for clarity
    kappaQ = param(1);
    thetaQ = param(2);
    eta = param(3);
    rho = param(4);
    V0 = param(5);

    flag = 0; %flag if the parameters are out of wanted intervals
    loss =0; %initialize loss

    if any(param>var_max) %check if any params are larger than they should be
        loss = loss+1e6; %Punish the loss heavily in that case
        flag = 1; % give flag out of pounds number

    end
    %same for the minimium values
    if any(param< var_min)
        loss = loss+1e6;
        flag = 1;

    end
    %Did I want to see the parameters when I call the function?
    if settings.displayParameters == 1
        disp([kappaQ, thetaQ, eta, rho, V0]);
    end
    
    %put the params back in an array in the wanted order
    parameters = {V0, thetaQ, kappaQ, eta, rho};
    
    %initialize matrices for call prices and volatilities
    callPrices_matrix = [];
    volatility_matrix = [];

    if flag == 0 %If parameters are in the acceptable range of values. If not dont waste time on these calculations

        for time = 1:length(T) % the functions dont take time as matrix for input so calculating for each time step

            callPrices = S0.*CallPricingFFT(model, n, 1, K./S0, T(time), r, 0, parameters{:}); %Calculate call prices using the specified model and paramters with CallPricingFFT.m
            callPrices(callPrices <= 0) = NaN; %if some values are 0 or less than zero make the NaN, since they dont make sense assuming no arbitrage
            impvol = blsimpv(S0,K,r,T(time),callPrices); %calculate the black-scholes implied volatilities from the prices

            %store call prices and implied volatilites
            callPrices_matrix(time,:) = callPrices; % not necessery, but usefull to check if something is odd
            volatility_matrix(time,:) = impvol;
        end



        volatility_matrix = fillmissing(volatility_matrix,'linear'); %interpolation of NaN values with matlab fillmissing values function
        %calculate the loss with immse. A matlab function to calculate mse of 2 matrix
       % loss = loss+immse(IVolSurf, volatility_matrix);
        %for some reason there is a difference between these 2 losses and I
        %get better results with immse
        square = (IVolSurf-volatility_matrix).^2;
        loss = loss+mean(square(:));
       %Error on the instructions is Squared error and not MSE, but in the
       %slides when talking about optimization all of the errors are MSE
       %So for this I assumed that IV-MSE should be used instead of the
       %squared error
       %I got better results with With IMMSE and TolFun = 1e-6;TolX =
       %1e-6;, but I dont know if thats "legal" to use on this project :/
        

        
        if settings.indPlotSurface == 1 % do I want to plot the volatility surfaces?
           
            figure(1); 
            surf(volatility_matrix); % Model output
            hold on;
            surf(IVolSurf);          % Target output
            title('Fitting Volatility Surface');
            xlabel('Strike Price');
            ylabel('Time to Maturity');
            zlabel('Implied Volatility');
            legend('Model output', 'Target output');
            hold off;
        end
        if settings.standardErrors == true %If I want to calculate standard errors using the jacobian.m it requires the model output of the volatility matrix
            loss = volatility_matrix(:);
        end

    end

end