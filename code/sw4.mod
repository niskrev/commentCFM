//=========================================================================
// SW(2007) AER model, 4 structural schoks (TFP, inv, gov't, MP) 
//=========================================================================

//=========================================================================
// ENDOGENOUS VARIABLES
var  
zcapf rkf kf pkf cf invef yf labf wf rrf mc zcap rk k pk 
c inve y lab pinf w r a g qs kpf kp rtil;    
// measurables
var
y_o c_o i_o w_o pi_o r_o h_o;
 

//=========================================================================
// EXOGENOUS INNOVATIONS
varexo 
ea eg eqs em;  
 
//=========================================================================
// DEEP AND STUCTURAL PARAMETERS
parameters 
curvw curvp cgy calfa czcap cbeta csadjcost ctou csigma chabb cfc cindw 
cprobw cindp cprobp csigl clandaw crpi crdy cry crr crhoa crhog 
crhoqs cg cgamma  
//clandap cbetabar crk cw cikbar cik clk cky ciy ccy crkky cwhlc cwly
std_a std_g std_qs std_m;
// steady states parameters
//y_ss c_ss i_ss w_ss pi_ss r_ss h_ss; 

//=========================================================================
// fixed parameters (these are the parameters that are fixed in SW(2007))
ctou    = 0.025;    // depreciation rate
curvp   = 10;      // good markets kimball aggragator 
curvw   = 10;      // labor markets kimball aggragator 
clandaw = 1.5;     // elasticity of substitution labor
cg      = 0.18;    // gov't consumption output share

//=========================================================================
// Deep Parameters values (these parameters are estimated in SW(2007))
// phip is one plus the share of fixed costs in production,
 
cbeta       = 0.99;     // time discount factor
cfc         = 1.61;     // (PHI) in the table 1A (AER 2007),
csadjcost   = 5.48;     // (phi) in the table 1A (AER 2007), is the steady-state elasticity of the capital adjustment cost function

calfa       = 0.20;      // (alpha) capital share
chabb       = 0.71;      // (h) Habit in consumtion     
cprobw      = 0.73;      // (xi_w) wage stickiness
cprobp      = 0.65;      // (xi_p) price stickiness
cindw       = 0.59;      // (i_w) wage indexation
cindp       = 0.47;      // (i_p) price indexation
csigl       = 1.92;      // (sigma_l) inverse of frish elasticity
csigma      = 1.39;      // (sigma_c) elasticity of substitution
czcap       = 0.26;      // ( czcap=(1-psi)/psi ) in the table 1A (AER 2007), is a positive
                         // function of the elasticity of the capital utilization 
                         // adjustment cost function and normalized to be between 0 1.
crpi        = 2.03;      // (r_{\pi}) taylor rules response to inflation
crr         = 0.87;      // (rho} )  taylor rules smoothing
cry         = 0.08;      // (r_{y}) taylor rules response to output gap
crdy        = 0.22;      // (r_{Dy}) taylor rules response to oputput growth

cgamma      = 1.000;     // trend, if =1 no trend

//=========================================================================
// Autoregressive paramters:
crhoa       = 0.95;      // (rho_a) autoregressive technology
crhog       = 0.97;      // (rho_g) autoregressive gov't
crhoqs      = 0.71;      // (rho_I) autoregressive investment
cgy         = 0.51;      // cross coefficient in the gov't process 

std_a  = 0.4618;
std_g  = 0.6090;
std_qs = 0.6017;
std_m  = 0.2513;

//=========================================================================
//=========================================================================
// SW MODEL LINEARIZED EQUATIONS
model(linear); 
//=========================================================================
// Composite paramters derived from steady state
# clandap     =  cfc;
# cbetabar    = cbeta*cgamma^(-csigma);
# crk         = (cbeta^(-1))*(cgamma^csigma) - (1-ctou);
# cw          = (calfa^calfa*(1-calfa)^(1-calfa)/(clandap*crk^calfa))^(1/(1-calfa));
# cikbar      = (1-(1-ctou)/cgamma);
# cik         = (1-(1-ctou)/cgamma)*cgamma;
# clk         = ((1-calfa)/calfa)*(crk/cw);
# cky         = cfc*(clk)^(calfa-1);
# ciy         = cik*cky;
# ccy         = 1-cg-cik*cky;
# crkky       = crk*cky;
# cwhlc       = (1/clandaw)*(1-calfa)/calfa*crk*cky/ccy;
# cwly        = 1-crk*cky;
//=========================================================================
// Steady states
# pi_ss  = 1;
# r_ss   = pi_ss/cbetabar;
# w_ss   = ( calfa^calfa * (1-calfa)^(1-calfa) / (1+clandaw) * crk^calfa )^( 1/(1-calfa) ) ;
# kh_ss  = calfa / (1-calfa) * w_ss / crk ;
# ky_ss  = csadjcost * (1/kh_ss)^(calfa-1);
# ik_ss  = (cgamma - (1-ctou));
# cy_ss  = 1 - cg - ik_ss * ky_ss;   
# ch_ss  = cy_ss / ky_ss * kh_ss ;
# h_ss   = (w_ss / (1+clandaw) * 1/(1-chabb/cgamma) * 1/ch_ss)^( 1/(csigl-1) ) ;
# k_ss   = kh_ss * h_ss;
# y_ss   = k_ss / ky_ss ;
# c_ss   = y_ss * cy_ss;
# i_ss   = k_ss * ik_ss;
//=========================================================================
// flexible economy
0*(1-calfa)*a + 1*a =  calfa*rkf+(1-calfa)*(wf)  ;
zcapf =  (1/(czcap/(1-czcap)))* rkf  ;
rkf   =  (wf)+labf-kf ;
kf    =  kpf(-1)+zcapf ;
invef = (1/(1+cbetabar*cgamma))* (  invef(-1) + cbetabar*cgamma*invef(1)+(1/(cgamma^2*csadjcost))*pkf ) +qs ;
pkf   = -rrf +(crk/(crk+(1-ctou)))*rkf(1) +  ((1-ctou)/(crk+(1-ctou)))*pkf(1) ;
cf    = (chabb/cgamma)/(1+chabb/cgamma)*cf(-1) + (1/(1+chabb/cgamma))*cf(+1) +((csigma-1)*cwhlc/(csigma*(1+chabb/cgamma)))*(labf-labf(+1)) - (1-chabb/cgamma)/(csigma*(1+chabb/cgamma))*(rrf) ;
yf    = ccy*cf+ciy*invef+g  +  crkky*zcapf ;
yf    = cfc*( calfa*kf+(1-calfa)*labf +a );
wf    = csigl*labf 	+(1/(1-chabb/cgamma))*cf - (chabb/cgamma)/(1-chabb/cgamma)*cf(-1) ;
kpf   = (1-cikbar)*kpf(-1)+(cikbar)*invef + (cikbar)*(cgamma^2*csadjcost)*qs ;

//=========================================================================
// sticky price - wage economy
mc   =  calfa*rk+(1-calfa)*(w) - 1*a - 0*(1-calfa)*a ;
zcap =  (1/(czcap/(1-czcap)))* rk ;
rk   = w + lab - k ;
k    = kp(-1)+zcap ;
kp   = (1-cikbar)*kp(-1)+cikbar*inve + cikbar*cgamma^2*csadjcost*qs ;
inve = (1/(1+cbetabar*cgamma))* (  inve(-1) + cbetabar*cgamma*inve(1) + 
       (1/(cgamma^2*csadjcost))*pk ) + qs ;
pk   = - r + pinf(1) +
       (crk/(crk+(1-ctou)))*rk(1) +  ((1-ctou)/(crk+(1-ctou)))*pk(1) ;
c    = (chabb/cgamma)/(1+chabb/cgamma)*c(-1) + (1/(1+chabb/cgamma))*c(+1) +
       ((csigma-1)*cwhlc/(csigma*(1+chabb/cgamma)))*(lab-lab(+1)) - 
       (1-chabb/cgamma)/(csigma*(1+chabb/cgamma))*(r-pinf(+1));
y    = ccy*c + ciy*inve + g + 1*crkky*zcap ;
y    = cfc*( calfa*k + (1-calfa)*lab + a );
pinf = (1/(1+cbetabar*cgamma*cindp)) * ( cbetabar*cgamma*pinf(1) + cindp*pinf(-1) + 
       ((1-cprobp)*(1-cbetabar*cgamma*cprobp)/ cprobp)/((cfc-1)*curvp+1)*(mc));
w    = (1/(1+cbetabar*cgamma))*w(-1) + 
       (cbetabar*cgamma/(1+cbetabar*cgamma))*w(1) + 
       (cindw/(1+cbetabar*cgamma))*pinf(-1) - 
       (1+cbetabar*cgamma*cindw)/(1+cbetabar*cgamma)*pinf + 
       (cbetabar*cgamma)/(1+cbetabar*cgamma)*pinf(1) + 
       (1-cprobw)*(1-cbetabar*cgamma*cprobw)/((1+cbetabar*cgamma)*cprobw) * 
       (1/((clandaw-1)*curvw+1)) * 
       (csigl*lab + (1/(1-chabb/cgamma))*c - ((chabb/cgamma)/(1-chabb/cgamma))*c(-1) -w);
r      = crpi*(1-crr)*pinf + cry*(1-crr)*(y-yf) + crdy*(y - yf) + crr*rtil(-1) + em*std_m;
rtil   = r - crdy/crr*(y-yf);
//=========================================================================
//  autoregressive processes
a  = crhoa*a(-1)  + ea*std_a;
g  = crhog*(g(-1)) + eg*std_g + cgy*ea*std_a;
qs = crhoqs*qs(-1) + eqs*std_qs;
 
//=========================================================================
//  observable variables:
y_o      = y + log(y_ss);
c_o      = c  + log(c_ss);
i_o      = inve + log(i_ss);
w_o      = w  + log(w_ss);
pi_o     = pinf + log(pi_ss);
r_o      = r + log(r_ss);
h_o      = lab + log(h_ss);

end; 
//=========================================================================
//=========================================================================

estimated_params;
ctou, 0.025;    // depreciation rate
//curvp, 10;      // good markets kimball aggragator 
//curvw, 10;      // labor markets kimball aggragator 
clandaw, 1.5;     // elasticity of substitution labor
cg, 0.18;      // gov't consumption output share
cbeta, 0.99;     // time discount factor
cfc, 1.61;     // (PHI) in the table 1A (AER 2007),
csadjcost, 5.48;     // (phi) in the table 1A (AER 2007), is the steady-state elasticity of the capital adjustment cost function
calfa, 0.20;      // (alpha) capital share
chabb, 0.71;      // (h) Habit in consumtion     
cprobw, 0.73;      // (xi_w) wage stickiness
cprobp, 0.65;      // (xi_p) price stickiness
cindw, 0.59;      // (i_w) wage indexation
cindp, 0.47;      // (i_p) price indexation
csigl, 1.92;      // (sigma_l) inverse of frish elasticity
csigma, 1.39;      // (sigma_c) elasticity of substitution
czcap, 0.26;      // ( czcap=(1-psi)/psi ) in the table 1A (AER 2007), is a positive
crpi , 2.03;      // (r_{\pi}) taylor rules response to inflation
crr  , 0.87;      // (rho} )  taylor rules smoothing
cry  , 0.08;      // (r_{y}) taylor rules response to output gap
crdy , 0.22;      // (r_{Dy}) taylor rules response to oputput growth
//cgamma, 1.000;     // trend, if =1 no trend
crhoa, 0.95;      // (rho_a) autoregressive technology
crhog, 0.97;      // (rho_g) autoregressive gov't
crhoqs, 0.71;      // (rho_I) autoregressive investment
cgy, 0.51;      // cross coefficient in the gov't process 
std_a, 0.4618;
std_g, 0.6090;
std_qs, 0.6017;
std_m, 0.2513;
end;

varobs h_o c_o i_o w_o; // pi_o r_o h_o;
//varobs y_o c_o i_o w_o pi_o r_o h_o;


shocks;
var ea;
stderr 0.4618;
var eg;
stderr 0.6090;
var eqs;
stderr 0.6017;
var em;
stderr 0.2513;
end;

//=========================================================================
check;
steady;
//stoch_simul(order=0,periods=200,noprint,nograph); 

//stoch_simul(irf=0,ar=10) y_o c_o i_o w_o pi_o r_o h_o;
 
identification;
