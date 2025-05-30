function settings = calibrationSettings
% In this function, the settings of the calibration problem are specified

%% GENERAL SETTINGS
% Number of days in a year
settings.tradingDaysInYear = 252;

%% CLOSED-FORM SOLUTION -SETTINGS
settings.n = 13;
settings.model = 'Heston';

%% INITIAL PARAMETER VALUES
kappa = 8;
theta = 0.25^2;
eta = 0.8;
rho = -0.8;
V0 = 0.3^2;

settings.parameters0 = [kappa; theta; eta; rho; V0];

%% MINUMUM AND MAXIMUM VALUES for the parameters of volatility model
settings.minKappa = 0.5; settings.maxKappa = 12;
settings.minTheta = 0.01^2; settings.maxTheta = 1;
settings.minEta = 0.05^2; settings.maxEta = 2;
settings.minRho = -1; settings.maxRho = 1;
settings.minV0 = 0.01^2; settings.maxV0 = 1;
settings.numberOfVariables = 5;

%% OPTIMIZATION SETTINGS
settings.calibrOptions.MaxFunEvals = 200*6;
settings.calibrOptions.MaxIter = 200*6;
settings.calibrOptions.TolFun = 1e-4;
settings.calibrOptions.TolX = 1e-4;
settings.calibrOptions.Display = 'iter';
settings.calibrOptions.FunValCheck = 'on';

%% DISPLAY SETTINGS, provisional result
settings.displayProvisionalResults = true;

if settings.displayProvisionalResults
    settings.indPlotSurface = 1;
    settings.displayParameters = 1;
    settings.calibrOptions.Display = 'iter';
    settings.calibrOptions.FunValCheck = 'on';
else
    settings.indPlotSurface = 0;
    settings.showParameters = 0;
    settings.calibrOptions.Display = 'off';
    settings.calibrOptions.FunValCheck = 'off';
end


end

