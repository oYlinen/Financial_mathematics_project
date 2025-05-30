% This is material illustrating the methods from the book
% Financial Modelling  - Theory, Implementation and Practice with Matlab
% source
% Wiley Finance Series
% ISBN 978-0-470-74489-5
%
% Date: 02.05.2012
%
% Authors:  Joerg Kienitz
%           Daniel Wetterau
%
% Please send comments, suggestions, bugs etc. to neueemail
%
% The authors are not responsible for any loss or damage from applying
% this code
%
% (C) Kienitz, Wetterau


function y = CharacteristicFunctionLib(model,u,lnS,T,r,d,varargin)
%---------------------------------------------------------
% Characteristic Function Library of the following models:
%---------------------------------------------------------
% Black Scholes
% Merton Jump Diffusion
% Heston Stochastic Volatility Model
% Bates Stochastic Volatility / Jump Diffusion Model
% Displaced Heston
% Heston-Hull-White
% Variance Gamma
% Normal Inverse Gaussian
% Meixner
% Generalized Hyperbolic
% CGMY
% Variance Gamma with Gamma Ornstein Uhlenbeck clock
% Variance Gamma with CIR clock
% NIG with Gamma Ornstein Uhlenbeck clock
% NIG with CIR clock
% SV with Levy jumps:
% SVJJ
% SVVG; stochastic volatility model with variance gamma Lévy jump in return
% SVLS: stochastic volatility model with log stable Lévy jump in return
% for further details about SVVG and SVLS, please refer to 'Return dynamics
% with Lévy jumps: evidence from stock and option prices' written by H.Li,
% et al.
% BNS model with Gamma SV(BNSSG)
% BNS model with IG SV(BNSSIG)
% Time-changed Levy models:
% NIGSA
% VGSA
% CGMYSA
% NIGSAM
% VGSAM
% CGMYSAM
% NIGSG
% VGSG
% NIGIG
% VGIG
% NIGSIG
% VGSIG
% for further details about time-changed model, please refer to 'Stochastic
% volatility for Lévy processes' written by P.Carr,
% et al.
%---------------------------------------------------------
%
optAlfaCalculation = true;

ME1 = MException('VerifyInput:InvalidNrOfArguments',...
    'Invalid number of Input arguments');
ME2 = MException('VerifyInput:InvalidModel',...
    'Undefined Model');

if strcmp(model,'BlackScholes')
    if nargin == 7
        funobj = @cf_bs;
    else 
        throw(ME1)
    end
elseif strcmp(model,'DDHeston')
    if nargin == 13
        funobj = @cf_ddheston;
    else 
        throw(ME1)
    end
elseif strcmp(model,'Heston')
    if nargin == 11
        funobj = @cf_heston;
    else 
        throw(ME1)
    end
elseif strcmp(model,'HestonHullWhite')
   if nargin == 14
       funobj = @cf_hestonhullwhite;
   else
       throw(ME1)
   end
elseif strcmp(model,'HHW')
   if nargin == 15
       funobj = @cf_hhw;
   else
       throw(ME1)
   end
elseif strcmp(model,'H1HW')
   if nargin == 15
       funobj = @cf_h1hw;
   else
       throw(ME1)
   end
elseif strcmp(model,'Merton')
    if nargin == 10
        funobj = @cf_merton;
    else 
        throw(ME1)
    end
elseif strcmp(model,'Bates')
    if nargin == 14
        funobj = @cf_bates;
    else 
        throw(ME1)
    end
elseif strcmp(model,'BatesHullWhite')
    if nargin == 17
        funobj = @cf_bateshullwhite;
    else 
        throw(ME1)
    end
elseif strcmp(model,'VarianceGamma')
    if nargin == 9
        funobj = @cf_vg;
    else 
        throw(ME1)
    end
elseif strcmp(model,'NIG')
    if nargin == 10
        funobj = @cf_nig;
    else
        throw(ME1)
    end
elseif strcmp(model,'CGMY')
    if nargin == 10
        funobj = @cf_cgmy;
    else
        throw(ME1)
    end
elseif strcmp(model,'Meixner')
    if nargin == 9
        funobj = @cf_meixner;
    else
        throw(ME1)
    end
elseif strcmp(model,'GH')
    if nargin == 10
        funobj = @cf_gh;
    else
        throw(ME1)
    end
elseif strcmp(model,'VarianceGammaCIR')
    if nargin == 12
        funobj = @cf_vg_cir;
    else
        throw(ME1)
    end
elseif strcmp(model,'VarianceGammaOU')
    if nargin == 12
        funobj = @cf_vg_gou;
    else
        throw(ME1)
    end
elseif strcmp(model,'NIGOU')
    if nargin == 12
        funobj = @cf_nig_gou;
    else
        throw(ME1)
    end
elseif strcmp(model,'NIGCIR')
    if nargin == 12
        funobj = @cf_nig_cir;
    else
        throw(ME1)
    end
elseif strcmp(model,'SVJJ')
    if nargin == 16
        funobj = @cf_svjj;
    else 
        throw(ME1)
    end
elseif strcmp(model,'SVJJ_P')
    if nargin == 18
        funobj = @cf_svjj_p;
    else 
        throw(ME1)
    end    
elseif strcmp(model,'SVVG')
    %if nargin == 15
    if nargin == 14
        funobj = @cf_svvg;
    else 
        throw(ME1)
    end
elseif strcmp(model,'SVNIG')
    %if nargin == 15
    if nargin == 14
        funobj = @cf_svnig;
    else 
        throw(ME1)
    end    
elseif strcmp(model,'SValphas')
    if nargin == 17
        funobj = @cf_svalphas;
    else 
        throw(ME1)
    end     
elseif strcmp(model,'BNSSG')
    if nargin == 11
        funobj = @cf_bnssg;
    else 
        throw(ME1)
    end
elseif strcmp(model,'BNSSIG')
    if nargin == 11
        funobj = @cf_bnssig;
    else 
        throw(ME1)
    end
elseif strcmp(model,'NIGSA')
    if nargin == 12
        funobj = @cf_nigsa;
    else 
        throw(ME1)
    end
elseif strcmp(model,'VGSA')
    if nargin == 12
        funobj = @cf_vgsa;
    else 
        throw(ME1)
    end
elseif strcmp(model,'CGMYSA')
    if nargin == 15
        funobj = @cf_cgmysa;
    else 
        throw(ME1)
    end
elseif strcmp(model,'NIGSAM')
    if nargin == 12
        funobj = @cf_nigsam;
    else 
        throw(ME1)
    end
elseif strcmp(model,'VGSAM')
    if nargin == 12
        funobj = @cf_vgsam;
    else 
        throw(ME1)
    end
elseif strcmp(model,'CGMYSAM')
    if nargin == 15
        funobj = @cf_cgmysam;
    else 
        throw(ME1)
    end
elseif strcmp(model,'NIGSG')
    if nargin == 14
        funobj = @cf_nigsg;
    else 
        throw(ME1)
    end
elseif strcmp(model,'VGSG')
    if nargin == 14
        funobj = @cf_vgsg;
    else 
        throw(ME1)
    end
elseif strcmp(model,'NIGIG')
    if nargin == 12
        funobj = @cf_nigig;
    else 
        throw(ME1)
    end
elseif strcmp(model,'VGIG')
    if nargin == 13
        funobj = @cf_vgig;
    else 
        throw(ME1)
    end
elseif strcmp(model,'NIGSIG')
    if nargin == 13
        funobj = @cf_nigsig;
    else 
        throw(ME1)
    end
elseif strcmp(model,'VGSIG')
    if nargin == 13
        funobj = @cf_vgsig;
    else 
        throw(ME1)
    end
else
    throw(ME2)
end

fval = feval(funobj,u,lnS,T,r,d,varargin{:});

if optAlfaCalculation == true
    y = fval;
else
    y = exp(fval);
end

end


%% Explicit Implementation of the characteristic Functions E[exp(iu*lnS_T)]
%-----------------------------------------------------------------------
   
function y = cf_bs(u,lnS,T,r,d,sigma)
% Black Scholes
    y = 1i*u*(lnS+(r-d-0.5*sigma*sigma)*T) - 0.5*sigma*sigma*u.*u*T;
end


function y = cf_merton(u,lnS,T,r,d,sigma,a,b,lambda)
% Merton Jump Diffusion
    y = cf_bs(u,lnS,T,r,d,sigma) ...
        + cf_jumplognormal(u,a,b,lambda,T);
end

function y = cf_ddheston(u,lnS,T,r,d,V0,theta,kappa,omega,rho,lambda,b)
% Displaced Diffusion Heston
U = 1i*u;
v = 0.5*(lambda*b)^2*U.*(U-1);
theta_star = kappa -rho*omega*lambda*b*U;
gamma = sqrt(theta_star.^2-2*omega^2*v);
Avu = kappa * theta / omega^2 *(2*log(2*gamma./(theta_star ...
    + gamma-exp(-gamma*T).*(theta_star-gamma)))+(theta_star-gamma)*T);
Bvu = 2*v.*(1-exp(-gamma*T))./((theta_star+gamma)...
    .*(1-exp(-gamma*T))+2*gamma.*exp(-gamma*T));

y = Avu + Bvu*V0 + U *(lnS+(r-d)*T);
end


function y = cf_heston(u,lnS,T,r,d,V0,theta,kappa,omega,rho)
% Heston  
alfa = -.5*(u.*u + u*1i);
beta = kappa - rho*omega*u*1i;
omega2 = omega * omega;
gamma = .5 * omega2;

D = sqrt(beta .* beta - 4.0 * alfa .* gamma);

bD = beta - D;
eDt = exp(- D * T);

G = bD ./ (beta + D);
B = (bD ./ omega2) .* ((1.0 - eDt) ./ (1.0 - G .* eDt));
psi = (G .* eDt - 1.0) ./(G - 1.0);
A = ((kappa * theta) / (omega2)) * (bD * T - 2.0 * log(psi));

y = A + B*V0 + 1i*u*(lnS+(r-d)*T);

end

function y = cf_hhw(u,lnS,T,r0,d,V0,theta,kappa,omega,lambda,eta,rho12,rho13,ircurve)
% Heston Hull White with % correlation(variance,rate) = 0
% dr(t) = lambda(r-r(t))dt + eta dW(t); r constant
    
    D1 = sqrt((omega*rho12*1i*u-kappa).^2-omega^2*1i*u.*(1i*u-1));
    g = (kappa-omega*rho12*1i*u-D1)./(kappa-omega*rho12*1i*u+D1);
    
    a = sqrt(theta - .125 * omega^2/kappa);
    b = sqrt(V0) - a;
    ct=.25*omega^2*(1-exp(-kappa))/kappa;
    lambdat=4*kappa*V0*exp(-kappa)/(omega^2*(1-exp(-kappa)));
    d2=4*kappa*theta/omega^2;
    F1 = sqrt(ct*(lambdat-1)+ct*d2+ct*d2/(2*(d2+lambdat)));
    c = -log((F1-a)/b);
    
    I2 = kappa*theta/omega^2*(T*(kappa-omega*rho12*1i*u-D1)-2*log((1-g.*exp(-D1*T))./(1-g)));
    I3 = eta^2*(1i+u).^2/(4*lambda^3)*(3+exp(-2*lambda*T)-4*exp(-lambda*T)-2*lambda*T);
    I4 = -eta*rho13/lambda *(1i*u+u.^2)*(b/c*(1-exp(-c*T))+a*T+a/lambda*(exp(-lambda*T)-1)+b/(c-lambda)*exp(-c*T)*(1-exp(-T*(lambda-c))));
    
    % curve stuff
    date_T = add2date(ircurve.Settle,T);
    Theta = (1-1i*u) * (log(ircurve.getDiscountFactors(date_T))+eta^2/(2*lambda^3)*(T/lambda+2*(exp(-lambda*T)-1)-0.5*(exp(-2*lambda)-1)));
    
    A = I2+I3+I4+Theta;
    BV = (1-exp(-D1*T))./(omega^2.*(1-g.*exp(-D1*T))).*(kappa-omega*rho12*1i*u-D1);
    Br = (1i*u-1)/lambda*(1-exp(-lambda*T));
    
    y = A + 1i*u * (lnS + (r0-d)*T) + BV * V0  + Br * r0;
end


function y = cf_hullwhite(u,T,lambda,eta,ircurve)
% Hull White  
    %maturity dates
    date_t = add2date(ircurve.Settle,0);
    date_T = add2date(ircurve.Settle,T);
 
    %time to maturity
    tau_0_T = diag(yearfrac(ircurve.Settle,date_T));
    %CURRENTLY equals zero
    tau_0_t = diag(yearfrac(ircurve.Settle,date_t));
    %CURRENTLY equal tau_0_T
    tau_t_T = diag(yearfrac(date_t,date_T));
    
    %used for short rate and  instantaneous forward rate calculations
    Delta = 1e-6;
    P0_T = ircurve.getDiscountFactors(date_T);
    %CURRENTLY P(0,t) equals one since t = 0
    P0_t = 1.0;
    P0_t_plus_Delta = interp1([ircurve.Settle;ircurve.Dates],...
        [1.0;ircurve.Data],datenum(date_t+Delta),'linear');

    %CURRENTLY short rate = instantaneous forward rate
    %--------------------------------------------------
    %short rate
    instR = (1.0/P0_t_plus_Delta-1)/yearfrac(ircurve.Settle,datenum(date_t+Delta));
    %instantaneous forward rate
    instF = -(P0_t_plus_Delta - P0_t)/yearfrac(ircurve.Settle,datenum(date_t+Delta))./P0_t;
    
    %CURRENTLY not used
    %----------------------------------------------
    %dt = yearfrac(ircurve.Settle, ircurve.Dates);
    %fwd_dates = add2date(date_t,dt);
    %fwd_rates = ircurve.getForwardRates(fwd_dates);
    %P0_t_T = ircurve.getDiscountFactors(date_T)./ircurve.getDiscountFactors(date_t);
    %P0_t_T_plus_Delta = ircurve.getDiscountFactors(datenum(date_T+Delta))./ircurve.getDiscountFactors(datenum(date_t));
    %------------------------------------------------

    %auxiliary variables
    aux = eta*eta/lambda/lambda;
    sigma2_func = @(t)aux*(t + (exp(-lambda*t)...
        .*(2.0 - 0.5*exp(-lambda*t))-1.5)/lambda);
    B_t_T = (1.0-exp(-lambda*tau_t_T))/lambda;
    %CURRENTLY psi_t = instantaneous forward rate since tau_0_t = 0
    psi_t = instF + .5*aux*(1-exp(-lambda*tau_0_t)).^2;

    %variance of integrated short rate
    sigma2_R = feval(sigma2_func,tau_t_T);
    
    %mean of integrated short rate
    mu_R = B_t_T.*(instR - psi_t) + diag(log(P0_t./P0_T)) ...
        + 0.5*(feval(sigma2_func,tau_0_T) - feval(sigma2_func,tau_0_t));
    
    y = 1i*u*mu_R - 0.5*u.*u*sigma2_R;
end

function y = cf_hestonhullwhite(u,lnS,T,r,d,V0,theta,kappa,omega,rho,lambda,eta,ircurve)
% Heston Hull White with correlation(variance,rate)=correlation(asset,rate)=0
% dr(t) = lambda(curve-r(t))dt + eta dW(t); curve is intial term structure
y = cf_heston(u,lnS,T,0,d,V0,theta,kappa,omega,rho) ...
        + cf_hullwhite(u+1i,T,lambda,eta,ircurve);
end

function y = cf_h1hw(u,lnS,T,r0,d,V0,theta,kappa,omega,lambda,eta,rho12,rho13,thetar)
% Heston Hull White with % correlation(variance,rate) = 0
% dr(t) = lambda(r-r(t))dt + eta dW(t); r constant
    
    D1 = sqrt((omega*rho12*1i*u-kappa).^2-omega^2*1i*u.*(1i*u-1));
    g = (kappa-omega*rho12*1i*u-D1)./(kappa-omega*rho12*1i*u+D1);
    
    a = sqrt(theta - .125 * omega^2/kappa);
    b = sqrt(V0) - a;
    ct=.25*omega^2*(1-exp(-kappa))/kappa;
    lambdat=4*kappa*V0*exp(-kappa)/(omega^2*(1-exp(-kappa)));
    d2=4*kappa*theta/omega^2;
    F1 = sqrt(ct*(lambdat-1)+ct*d2+ct*d2/(2*(d2+lambdat)));
    c = -log((F1-a)/b);
    
    I1 = thetar * (1i*u-1) * (T+(exp(-lambda * T)-1)/lambda);
    I2 = kappa*theta/omega^2*(T*(kappa-omega*rho12*1i*u-D1)-2*log((1-g.*exp(-D1*T))./(1-g)));
    I3 = eta^2*(1i+u).^2/(4*lambda^3)*(3+exp(-2*lambda*T)-4*exp(-lambda*T)-2*lambda*T);
    I4 = -eta*rho13/lambda *(1i*u+u.^2)*(b/c*(1-exp(-c*T))+a*T+a/lambda*(exp(-lambda*T)-1)+b/(c-lambda)*exp(-c*T)*(1-exp(-T*(lambda-c))));

    
    A = I1+I2+I3+I4;
    BV = (1-exp(-D1*T))./(omega^2.*(1-g.*exp(-D1*T))).*(kappa-omega*rho12*1i*u-D1);
    Br = (1i*u-1)/lambda*(1-exp(-lambda*T));
    
    y = A + 1i*u * (lnS + (r0-d)*T) + BV * V0  + Br * r0;
end

function y = cf_bates(u,lnS,T,r,d,V0,theta,kappa,omega,rho,a,b,lambda)
% Bates
    y = cf_heston(u,lnS,T,r,d,V0,theta,kappa,omega,rho)...
    + cf_jumplognormal(u,a,b,lambda,T);
end

function y = cf_bateshullwhite(u,lnS,T,r,d,V0,theta,kappa,omega,rho,a,b,lambda,lambda1,eta,ircurve)
% Bates Hull White
    y = cf_heston(u,lnS,T,r,d,V0,theta,kappa,omega,rho)...
        + cf_jumplognormal(u,a,b,lambda,T)...
        + cf_hullwhite(u+1i,T,lambda1,eta,ircurve);
end

function yJump = cf_jumplognormal(u,a,b,lambda,T)
% LogNormalJump for Merton and Bates
    yJump = lambda*T*(-a*u*1i + (exp(u*1i*log(1.0+a)+0.5*b*b*u*1i.*(u*1i-1.0))-1.0));
end

function y = cf_vg(u,lnS,T,r,d,sigma,nu,theta)
% Variance Gamma
    omega = (1/nu)*( log(1-theta*nu-sigma*sigma*nu/2) );
    tmp = 1 - 1i * theta * nu * u + 0.5 * sigma * sigma * u .* u * nu;
    y = 1i * u * (lnS + (r + omega - d) * T ) - T*log(tmp)/nu;
end


function y = cf_nig(u,lnS,T,r,d,alfa,beta,mu,delta)
% Normal Inverse Gaussian
    m = delta*(sqrt(alfa*alfa-(beta+1)^2)-sqrt(alfa*alfa-beta*beta));
    tmp = 1i*u*mu*T-delta*T*(sqrt(alfa*alfa-(beta+1i*u).^2)...
        -sqrt(alfa*alfa-beta*beta));
    y = 1i*u*(lnS + (r-d+m)*T) + tmp;
end


function y = cf_meixner(u,lnS,T,r,d,alfa,beta,delta)
% Meixner
    m = -2*delta*(log(cos(0.5*beta)) - log(cos((alfa+beta)/2)));
    tmp = (cos(0.5*beta)./cosh(0.5*(alfa*u-1i*beta))).^(2*T*delta);
    y = 1i*u*(lnS + (r-d+m)*T) + log(tmp);
end


function y = cf_gh(u,lnS,T,r,d,alfa,beta,delta,nu)
% Generalized Hyperbolic
    arg1 = alfa*alfa-beta*beta;
    arg2 = arg1-2*1i*u*beta + u.*u;
    argm = arg1-2*beta-1;
    m = -log((arg1./argm).^(0.5*nu).*besselk(nu,delta*sqrt(argm))./besselk(nu,delta*sqrt(arg1)));
    tmp = (arg1./arg2).^(0.5*nu).*besselk(nu,delta*sqrt(arg2))./besselk(nu,delta*sqrt(arg1));
    y = 1i*u*(lnS + (r-d+m)*T) + log(tmp).*T;
end


function y = cf_cgmy(u,lnS,T,r,d,C,G,M,Y)
% CGMY
    m = -C*gamma(-Y)*((M-1)^Y-M^Y+(G+1)^Y-G^Y);
    tmp = C*T*gamma(-Y)*((M-1i*u).^Y-M^Y+(G+1i*u).^Y-G^Y);
    y = 1i*u*(lnS + (r-d+m)*T) + tmp;
end

function y = cf_vg_cir(u,lnS,T,r,d,C,G,M,kappa,eta,lambda)
% VG CIR
y0 = 1;
psiX_u = (-1i)*C*(log(G*M)-log(G*M+(M-G)*1i*u + u.*u));
psiX_i = (-1i)*C*(log(G*M)-log(G*M+(M-G) -1 ));

gamma_u = sqrt(kappa^2-2*lambda^2*1i*psiX_u);
gamma_i = sqrt(kappa^2-2*lambda^2*1i*psiX_i); 

tmp = kappa^2*eta*T/lambda/lambda + 2*y0*1i*psiX_u./(kappa+gamma_u.*coth(gamma_u*T/2)) ...
    -log(cosh(0.5*gamma_u*T) + kappa*sinh(0.5*gamma_u*T)./gamma_u) ...
    *(2*kappa*eta/lambda/lambda);

y_T = kappa^2*eta*T/lambda/lambda + 2*y0*1i*psiX_i./(kappa+gamma_i*coth(gamma_i*T/2)) ...
     -log(cosh(0.5*gamma_i*T) + kappa*sinh(0.5*gamma_i*T)/gamma_i )...
        *(2*kappa*eta/lambda/lambda);

y = 1i*u*(lnS+(r-d).*T - y_T) + tmp;
end

function y = cf_vg_gou(u,lnS,T,r,d,C,G,M,lambda,a,b)
% VG GOU
y0 = 1;

psiX_u = (-1i)*C*(log(G*M)-log(G*M+(M-G)*1i*u + u.*u));
psiX_i = (-1i)*C*(log(G*M)-log(G*M+(M-G) -1 ));    
    
f2_u = 1i*psiX_u*T*lambda*a./(lambda*b-1i*psiX_u) ...
    + a*b*lambda./(b*lambda-1i*psiX_u)...
    .*log(1 - 1i*psiX_u/(lambda*b)*(1-exp(-T*lambda)));
f3_u = 1/lambda*psiX_u*(1-exp(-lambda*T));

f1_u = f2_u + 1i*y0*f3_u + a*log((1-1i/b*f3_u)./(1-1i/b*f3_u));

y_T = 1i*psiX_i*y0*lambda^(-1)*(1-exp(-lambda*T)) ...
    + lambda*a./(1i*psiX_i-lambda*b)...
    .*(b*log(b./(b-1i*psiX_i*lambda^(-1)*(1-exp(-lambda*T))))-1i*psiX_i*T);

y = 1i*u*(lnS + (r-d).*T - y_T) + f1_u;
end

function y = cf_nig_gou(u,lnS,T,r,d,alpha,beta,delta,lambda,a,b)
% NIG GOU
y0 = 1;

psiX_u = (-1i) * (-delta)*(sqrt(alpha^2-(beta+1i*u).^2)-sqrt(alpha^2-beta^2));
psiX_i = (-1i) * (-delta)*(sqrt(alpha^2-(beta+1)^2)-sqrt(alpha^2-beta^2));

f2_u = 1i*psiX_u*T*lambda*a./(lambda*b-1i*psiX_u) ...
    + a*b*lambda./(b*lambda-1i*psiX_u)...
    .*log(1 - 1i*psiX_u/(lambda*b)*(1-exp(-T*lambda)));
f3_u = 1/lambda*psiX_u*(1-exp(-lambda*T));

f1_u = f2_u + 1i*y0*f3_u + a*log((1-1i/b*f3_u)./(1-1i/b*f3_u));

y_T = 1i*psiX_i*y0*lambda^(-1)*(1-exp(-lambda*T)) ...
    + lambda*a./(1i*psiX_i-lambda*b)...
    .*(b*log(b./(b-1i*psiX_i*lambda^(-1)*(1-exp(-lambda*T))))-1i*psiX_i*T);

y = 1i*u*(lnS + (r-d).*T - y_T) + f1_u;

end


function y = cf_nig_cir(u,lnS,T,r,d,alpha,beta,delta,kappa,eta,lambda)
% NIG CIR
y0 = 1;

psiX_u = (-1i) * (-delta)*(sqrt(alpha^2-(beta+1i*u).^2)-sqrt(alpha^2-beta^2));
psiX_i = (-1i) * (-delta)*(sqrt(alpha^2-(beta+1)^2)-sqrt(alpha^2-beta^2));

gamma_u = sqrt(kappa^2-2*lambda^2*1i*psiX_u);
gamma_i = sqrt(kappa^2-2*lambda^2*1i*psiX_i); 

tmp = kappa^2*eta*T/lambda/lambda + 2*y0*1i*psiX_u./(kappa+gamma_u.*coth(gamma_u*T/2)) ...
    -log(cosh(0.5*gamma_u*T) + kappa*sinh(0.5*gamma_u*T)./gamma_u) ...
    *(2*kappa*eta/lambda/lambda);

y_T = kappa^2*eta*T/lambda/lambda + 2*y0*1i*psiX_i./(kappa+gamma_i*coth(gamma_i*T/2)) ...
     -log(cosh(0.5*gamma_i*T) + kappa*sinh(0.5*gamma_i*T)/gamma_i )...
        *(2*kappa*eta/lambda/lambda);

y = 1i*u*(lnS+(r-d).*T - y_T) + tmp;
end


function y=cf_svjj(v,lnS,T,r,d,V0,lambda,muQy,sigmay,rhoJ,muv,sigmav,kq,rho,thetaq)
% SVJJ: the physical drift mu=r-0.5*v_t+eta_s*v_t+ksi(-1i,lambda,muQy,sigmay,rhoJ,muv), where
% km=kq-1i.*u.*rho.*sigmav;
% e=1i.*u+u.^2;
% s=sqrt(km.^2+e.*(sigmav.^2));
% b=e.*(1-exp(-T.*s))./(s+km+(s-km).*exp(-T.*s));
% KSI=ksi_svjj(u,lambda,muQy,sigmay,rhoJ,muv);
% c=lambda.*T-kq.*thetaq.*T.*e./(s-km)-lambda.*T.*(s-km)./(s-km-muv.*e)+2*kq.*thetaq.*e.*log(-km+exp(s.*t).*km+s+exp(s.*t).*s)./(s.^2-km.^2)...
%     +2*lambda.*muv.*e.*log(-km+exp(s.*t).*km+s+s.*exp(s.*t)-muv.*e+exp(s.*t).*e.*muv)./(s.^2-km.^2-2*km.*muv.*e-muv.^2.*e.^2);
% y=(-c-b.*V0+1i*u.*lnS+1i*u.*(r-d+ksi_svjj(-1i,lambda,muQy,sigmay,rhoJ,muv)).*T)+(-KSI.*T);
BETA = kq - rho*sigmav*1i*v;
GAMMA = (sigmav^2)/2;
ALPHA = -0.5*(v.^2 + 1i*v);
e = sqrt(BETA.^2 - 4*ALPHA*GAMMA);
rneg = (BETA - e)./(sigmav^2);
rpos = (BETA + e)/(sigmav^2);
g = rneg./rpos;
C = kq * (rneg*T - (2/(sigmav^2)) * log((1 - g.*exp(-e*T))./ (1 - g)));
D = rneg .* ((1 - exp(-e*T)) ./ (1 - g.*exp(-e*T)));
MUJ = exp(muQy + 0.5*sigmay^2) / (1 - rhoJ*muv) - 1;
c = 1 - rhoJ*muv*1i*v;
nu = ( (BETA + e) ./ ((BETA + e).*c - 2*muv*ALPHA) ) * T + ...
( (4*muv*ALPHA) ./ ((e.*c).^2 - (2*muv*ALPHA - BETA.*c) ...
.^2) ) .* log( 1 - ( ((e-BETA).*c + 2*muv*ALPHA) ./ ...
(2*e.*c) ).*(1 - exp(-e*T)) );
P = -T*(1 + MUJ*1i*v) + exp( muQy*1i*v + 0.5*(sigmay^2)*(1i*v).^2 ).*nu;
y = C*thetaq + D*V0 + P*lambda + 1i*v*(lnS + (r-d)*T);
end

function y=cf_svjj_p(u,lnS,T,r,d,V0,lambda,mu,sigmay,rhoJ,muv,sigmav,k,etav,etas,rho,theta)
% lnS: log stock price
% r: interest rate
% d: dividend
% V0: variance 
% physical drift mu=r-0.5*v_t+eta_s*v_t+ksi(-1i,lambda,muQy,sigmay,rhoJ,muv)
% etav, etas are risk premia.
% lambda, sigmay, rhoJ, muv, sigmav, k, rho, and theta are parameters under
% physical measure. 
muQy=log((1-(mu-r+d+0.5*V0-etas*V0)/lambda)*(1-muv*rhoJ))-0.5*sigmay^2;
muQy = 0;
% muQy is the mean jump size for the jump in return under risk neutral
% measure
kq=k-etav; % speed of mean reversion under risk neutral measure
thetaq=k.*theta./kq;  % long-run mean of variance under risk neutral measure
KSI=ksi_svjj(u,lambda,muQy,sigmay,rhoJ,muv);
km=k-1i*sigmav.*rho.*u;
e=1i*u+u.^2;
s=sqrt(km.^2+(sigmav.^2).*e);
b=(e.*(1-exp(-T.*s)))./(s+km+(s-km).*exp(-T.*s));
J=(T.*(km-s).*(k+s+muv.*e)+2.*muv.*e.*log(-km+s-muv.*e+(km+s+muv.*e).*exp(T.*s)))./((km-s+muv.*e).*(km+s+muv.*e));
c=kq.*thetaq.*(2*log((2*s-(s-km).*(1-exp(-s.*T)))./(2*s))+(s-km).*T)./(sigmav.^2)+lambda.*T-lambda.*J;
y=(-c-b.*V0+1i*lnS.*u+1i*u.*(r-d+ksi_svjj(-1i,lambda,muQy,sigmay,rhoJ,muv)).*T)+(-KSI.*T);
end
function ksi=ksi_svjj(x,lambda,muQ_y,sigmay,rhoJ,muv)
ksi=lambda*(1-exp(1i*x.*muQ_y-(x.^2)*(sigmay.^2)/2)/(1-1i*x.*muv.*rhoJ));
end


function y=cf_svvg(u,lnS,T,r,d,V0,gamma,nu,sigma,kq,thetaq,sigmav,rho)
% stochastic volatility model with Variance gamma process as the jump in
% the return
KSI=ksi_svvg(u,gamma,nu,sigma);
km=kq-1i*u.*rho.*sigmav;
s=sqrt(km.^2+(1i*u+u.^2)*(sigmav.^2));
b=((1i*u+u.*u).*(1-exp(-T.*s)))./(s+km+(s-km).*exp(-T.*s));
c=kq.*thetaq.*(2*log((2*s-(s-km).*(1-exp(-s.*T)))./(2*s))+(s-km).*T)./(sigmav.^2);
y=(-c-b.*V0+1i*u.*lnS+1i*u.*(r-d+ksi_svvg(-1i,gamma,nu,sigma)).*T)+(-KSI.*T);
end
function ksi=ksi_svvg(x,gamma,nu,sigma)
temp=max(eps,1-1i*x.*gamma.*nu+0.5*(sigma.^2).*nu.*(x.^2));
ksi=log(temp)./nu;
end

function y=cf_svnig(u,lnS,T,r,d,V0,gamma,nu,sigma,kq,thetaq,sigmav,rho)
% stochastic volatility model with Variance gamma process as the jump in
% the return
KSI=ksi_svnig(u,gamma,nu,sigma);
km=kq-1i*u.*rho.*sigmav;
s=sqrt(km.^2+(1i*u+u.^2)*(sigmav.^2));
b=((1i*u+u.*u).*(1-exp(-T.*s)))./(s+km+(s-km).*exp(-T.*s));
c=kq.*thetaq.*(2*log((2*s-(s-km).*(1-exp(-s.*T)))./(2*s))+(s-km).*T)./(sigmav.^2);
y=(-c-b.*V0+1i*u.*lnS+1i*u.*(r-d+ksi_svnig(-1i,gamma,nu,sigma)).*T)+(-KSI.*T);
end
function ksi=ksi_svnig(x,gamma,nu,sigma)
if sigma~=0
temp=max(eps,nu.^2./(sigma.^2)+x.^2-2*1i*gamma.*x./(sigma.^2));
ksi=-nu+sigma.*sqrt(temp);
else
ksi=0;
end
end


function y=cf_svalphas(u,lnS,T,r,d,V0,alpha,beta,gamma,sigma,k,theta,sigmav,rho,etav)
% stochastic volatility model with alpha-stbale Levy process as the jump in
% the return. If beta=-1, gamma=0, the process reduces to log-stable
% process
KSI=ksi_svalphas(u,alpha,beta,gamma,sigma);
km=k-etav-1i*u.*rho.*sigmav;
s=sqrt(km.^2+(1i*u+u.^2)*(sigmav.^2));
b=((1i*u+u.*u).*(1-exp(-T.*s)))./(s+km+(s-km).*exp(-T.*s));
c=k.*theta.*(2*log((2*s-(s-km).*(1-exp(-s.*T)))./(2*s))+(s-km).*T)./(sigmav.^2);
y=(-c-b.*V0+1i*u.*lnS+1i*u.*(r-d+ksi_svalphas(-1i,alpha,beta,gamma,sigma)).*T)+(-KSI.*T);
end
function ksi=ksi_svalphas(x,alpha,beta,gamma,sigma)
ksi=(sigma.^alpha).*(abs(x).^alpha).*(1-1i*beta.*sign(x).*tan(0.5*pi*alpha))+1i*gamma.*x;
end


function y=cf_bnssg(u,lnS,T,r,d,V0,rho,lambda,a,b)
% BNS model with gamma SV
f1=1i*u.*rho-0.5*(u.*u+1i*u).*(1-exp(-lambda.*T));
f2=1i*u.*rho-0.5*(u.*u+1i*u);
y=1i*u.*(lnS+T.*(r-d-a*lambda.*rho./(b-rho)))-0.5*V0.*(u.^2+1i*u).*(1-exp(-lambda.*T))./lambda+a.*(f2.*lambda.*T+b.*log((b-f1)./(b-1i*u.*rho)))/(b-f2);
end

function y=cf_bnssig(u,lnS,T,r,d,V0,rho,lambda,a,b)
% BNS model with IG SV
f1=1i*u.*rho-0.5*(u.*u+1i*u).*(1-exp(-lambda.*T));
f2=1i*u.*rho-0.5*(u.*u+1i*u);
y=1i*u.*(lnS+T.*(r-d-rho.*lambda.a./b./sqrt(1-2*rho./b./b)))-0.5*V0.*(u.^2+1i*u).*(1-exp(-lambda.*T))./lambda+...
    a.*(sqrt(b.*b-2*f1)-sqrt(b.*b-2*1i*u.*rho))+2*a.*f2.*(arctan(sqrt((b.*b-2*1i*u.*rho)./(2*f2-b.*b)))-arctan(sqrt((b.*b-2*f1)./(2*f2-b.*b))))./sqrt(2*f2-b.*b);
end


function y=cf_nigsa(u,lnS,T,r,d,sigma,nu,theta,k,eta,lambda)
% NIGSA
y=1i*u.*(lnS+(r-d).*T)+log(phi_CIR(-1i*phi_NIG(u,1,nu,theta),T,sigma,k,eta,lambda)./(phi_CIR(-1i*phi_NIG(-1i,1,nu,theta),T,sigma,k,eta,lambda).^{1i*u}));
end
function phi=phi_NIG(u,sigma,nu,theta)
phi=sigma.*(nu./sigma-sqrt(nu.*nu./sigma./sigma-2*theta.*1i*u./sigma./sigma+u.*u));
end
function phi=phi_CIR(u,T,y,k,eta,lambda)
gamma=sqrt(k.*k-2*lambda.*lambda.*1i*u);
phi=exp(k.*k.*eta.*T./lambda./lambda+2*1i*u.*y./(k+gamma.*coth(gamma.*T./2)))./((cosh(gamma.*T./2)+k.*sinh(gamma.*T./2)/gamma).^(2*k.*eta./lambda./lambda));
end


function y=cf_vgsa(u,lnS,T,r,d,C,G,M,k,eta,lambda)
% VGSA
y=1i*u.*(lnS+(r-d).*T)+log(phi_CIR(-1i*phi_VG(u,1,G,M),T,C,k,eta,lambda)./(phi_CIR(-1i*phi_VG(-1i,1,G,M),T,C,k,eta,lambda)).^(1i*u));
end
function phi=phi_VG(u,C,G,M)
phi=C.*log(G.*M./(G.*M+(M-G).*1i*u+u.*u));
end


function y=cf_cgmysa(u,lnS,T,r,d,C,G,M,Yp,Yn,zeta,k,eta,lambda)
% CGMYSA
y=1i*u.*(lnS+(r-d).*T)+log(phi_CIR(-1i*phi_cgmy(u,1,G,M,Yp,Yn,zeta),T,C,k,eta,lambda)./((phi_CIR(-1i*phi_cgmy(-1i,1,G,M,Yp,Yn,zeta),T,C,k,eta,lambda)).^(1i*u)));
end
function phi=phi_cgmy(u,C,G,M,Yp,Yn,zeta)
phi=C.*(gamma(-Yp).*((M-1i*u).^Yp-M.^Yp)+zeta.*gamma(-Yn).*((G+1i*u).^Yn-G.^Yn));
end


function y=cf_nigsam(u,lnS,T,r,d,sigma,nu,theta,k,eta,lambda)
% NIGSAM
y=1i*u.*(lnS+(r-d).*T)+log(phi_CIR(-1i*phi_NIG(u,1,nu,theta)-u.*phi_NIG(-1i,1,nu,theta),T,sigma,k,eta,lambda));
end


function y=cf_vgsam(u,lnS,T,r,d,C,G,M,k,eta,lambda)
% VGSAM
y=1i*u.*(lnS+(r-d).*T)+log(phi_CIR(-1i*phi_VG(u,1,G,M)-u.*phi_VG(-1i,1,G,M),T,C,k,eta,lambda));
end


function y=cf_cgmysam(u,lnS,T,r,d,C,G,M,k,eta,lambda)
% CGMYSAM
y=1i*u.*(lnS+(r-d).*T)+log(phi_CIR(-1i*phi_cgmy(u,1,G,M,Yp,Yn,zeta)-u.*phi_cgmy(-1i,1,G,M,Yp,Yn,zeta),T,C,k,eta,lambda));
end


function y=cf_nigsg(u,lnS,T,r,d,sigma,nu,theta,k,lambda,zeta,rho)
% NIGSG
y=1i*u.*(lnS+(r-d).*T)+log(phi_SG(-1i*phi_NIG(u,sigma,nu,theta),rho.*u,T,k,lambda,zeta,sigma))-1i*u.*log(phi_SG(-1i*phi_NIG(-1i,sigma,nu,theta),-1i*rho,T,k,lambda,zeta,sigma));
end
function phi=phi_SG(a,b,T,k,lambda,zeta,y0)
phi=exp(1i*a.*y0.*(1-exp(-k.*T))./k).*exp(ksi_SG(b+a.*(1-exp(-k.*T))./k,a,b,k,lambda,zeta)-ksi_SG(b,a,b,k,lambda,zeta));
end
function ksi=ksi_SG(x,a,b,k,lambda,zeta)
ksi=log((x+1i/zeta).^(lambda./(k-1i.*zeta.*(a+k.*b))).*((a+k.*b-k.*x).^(lambda.*zeta.*(a+k.*b))./k./((a+k.*b).*zeta+1i*k)));
end


function y=cf_vgsg(u,lnS,T,r,d,C,G,M,k,lambda,zeta,rho)
% VGSG
y=1i*u.*(lnS+(r-d).*T)+log(phi_SG(-1i*phi_VG(u,C,G,M),rho.*u,T,k,lambda,zeta,C))-1i*u.*log(phi_SG(-1i*phi_VG(-1i,C,G,M),-1i*rho,T,k,lambda,zeta,C));
end


function y=cf_nigig(u,lnS,T,r,d,sigma,nu,theta,mu,rho)
% NIGIG
y=1i*u.*(lnS+(r-d).*T)+log(phi_IG(-1i*ksi_NIG(u,sigma,nu,theta),rho.*u,T,k,mu,C))-1i*u.*log(phi_IG(-1i*ksi_NIG(-1i,sigma,nu,theta),-1i*rho,T,k,mu,C));
end
function phi=phi_IG(a,b,T,k,mu,y0)
phi=exp(1i*a.*y0.*(1-exp(-k.*T))./k).*exp(ksi_IG(b+a.*(1-exp(-k.*T))./k,a,b,k,mu)-ksi_IG(b,a,b,k,mu));
end
function ksi=ksi_IG(x,a,b,k,mu)
ksi=2*sqrt(mu.*mu-2*1i*x)./k-mu.*log(a+b.*k-k.*x)./k+2*sqrt(mu.*mu.*k-2*1i*(a+k.*b)).*arctanh(sqrt(k.*(mu.*mu-2*1i*x))./sqrt(mu.*mu.*k-2*1i*(a+k.*b)))./(k.^1.5);
end


function y=cf_vgig(u,lnS,T,r,d,C,G,M,k,mu,rho)
% VGIG
y=1i*u.*(lnS+(r-d).*T)+log(phi_IG(-1i*ksi_VG(u,C,G,M),rho.*u,T,k,mu,C))-1i*u.*log(phi_IG(-1i*ksi_VG(-1i,C,G,M),-1i*rho,T,k,mu,C));
end


function y=cf_nigsig(u,lnS,T,r,d,sigma,nu,theta,k,mu,rho)
% NIGSIG
y=1i*u.*(lnS+(r-d).*T)+log(phi_SIG(-1i*ksi_NIG(u,sigma,nu,theta),rho.*u,T,k,mu,sigma))-1i*u.*log(phi_SIG(-1i*ksi_NIG(-1i,sigma,nu,theta),-1i*rho,T,k,mu,sigma));
end
function phi=phi_SIG(a,b,T,k,mu,y0)
phi=exp(1i*a.*y0.*(1-exp(-k.*T))./k).*exp(ksi_SIG(b+a.*(1-exp(-k.*T))./k,a,b,k,mu)-ksi_SIG(b,a,b,k,mu));
end
function ksi=ksi_SIG(x,a,b,k,mu)
ksi=sqrt(mu.*mu-2*1i*x)./k+2*1i*(a+k.*b)*arctanh(sqrt(k.*(mu.*mu-2*1i*x))./sqrt(mu.*mu.*k-2*1i*(a+k.*b)))./(k.^1.5)./sqrt(mu.*mu.*k-2*1i*(a+k.*b));
end


function y=cf_vgsig(u,lnS,T,r,d,C,G,M,k,mu,rho)
% VGSIG
y=1i*u.*(lnS+(r-d).*T)+log(phi_SIG(-1i*ksi_VG(u,C,G,M),rho.*u,T,k,mu,C))-1i*u.*log(phi_SIG(-1i*ksi_VG(-1i,C,G,M),-1i*rho,T,k,mu,C));
end
























 
