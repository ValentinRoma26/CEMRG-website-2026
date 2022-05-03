function land_rat_tako
 set(0,'defaulttextfontsize',15);
 set(0,'defaultaxesfontsize',15);  

    
 ca=[0.207 0.199 0.190 0.214 0.222 0.300 0.375 0.413 0.417 0.433 0.449 0.446 0.401 0.418 0.380 0.387 0.385 0.372 0.355 0.343 0.342 0.321 0.314 0.302 0.291 0.271 0.256 0.267 0.245 0.249 0.241 0.234 0.207 0.207 0.207];
 cabeta = [0.160 0.216 0.296 0.397 0.477 0.640 0.671 0.782 0.803 0.825 0.893 0.810 0.877 0.849 0.851 0.783 0.821 0.763 0.718 0.696 0.618 0.579 0.522 0.463 0.440 0.370 0.340 0.292 0.218 0.212 0.179 0.163 0.149 0.160 0.160 0.160];
 
 figure(1); plot(5*(0:length(ca)-1),ca,5*(0:length(cabeta)-1),cabeta,'LineWidth',3); xlim([0 167]);
 
 
 [t,y] = ode15s(@(t,y) landrat(t,y,ca,0.6),0:2000,landrat());
 [t,yb] = ode15s(@(t,y) landrat(t,y,cabeta,1),0:2000,landrat());
 
 figure(2); plot(t-2000+167,120*y(:,1),t-2000+167,120*yb(:,1),'LineWidth',3); xlim([0 167])
 
 figure(1); set(gca,'ytick',0:0.2:1);set(gca,'xtick',[0 167]); xlabel('Time (ms)'); ylabel('Ca^{2+} (\mu{}M)');
 legend('Control','\beta-stimulated   .','Location','Best');
 figure(2); set(gca,'ytick',0:20:80);set(gca,'xtick',[0 167]); xlabel('Time (ms)'); ylabel('Active tension (kPa)');

       
 for q=1:2
   
 for i=1:length(ca)-1
   b=(ca(i+1)-ca(i))/5;
   a = ca(i); %-b*5*(i-1); % a+b*5*(i-1) == ca(i) ->   a == ca(i)-b*5*(i-1)
   fprintf('  case mtime <= %d{dimensionless}: %.3f{uM}+(%.3f{uM})*(mtime-%d{dimensionless});\n',5*i,a,b,5*(i-1))
 end
  
 fprintf('  otherwise: %.3f{uM};\n\n\n',ca(1));
 ca=cabeta;
 end
 
 


function [dydt, Tension] = landrat(t,y,Cai,ca50p)
if nargin==0 || isempty(t)  
   dydt=[2e-03 8e-02 0 0]; % PACED NOKO
  return
end

dlambda_dt = 0;
lambda = 1;

Cai = interp1(5*(0:length(Cai)-1),Cai, mod(t,167) );

dydt = zeros(size(y));
%-------------------------------------------------------------------------------
% Parameters
%-------------------------------------------------------------------------------


if nargin >= 4
  ca50 = ca50p;
else
  ca50 = 0.7;    % MODEL DEFAULT = 0.8, SET TO CA
end

Tref = 120;
perm50 = 0.35;
TRPN_n = 2;
koff    = 0.1;
nperm = 5;
kxb = 0.02;
beta_ca = -1.5;
beta_0 = 1.65;


% State Variables
XB   = y(1);
TRPN = max(0,y(2));
Q1 = y(3);
Q2 = y(4);
Q = Q1+Q2;

%-------------------------------------------------------------------------------
% XB model
%-------------------------------------------------------------------------------


lambda = min(1.2,lambda);


ca50= ca50 * (1 + beta_ca*(lambda-1));
dTRPN  = koff * ( (Cai/ca50)^TRPN_n * (1-TRPN) - TRPN);


permtot = sqrt((TRPN/perm50)^nperm);
inprmt= min(100, 1/permtot );

dXB  = kxb*(permtot*(1-XB) - inprmt*XB);


dydt(1)  = dXB;
dydt(2)  = dTRPN;


%-------------------------------------------------------------------------------
% Velocity dependence
%-------------------------------------------------------------------------------

a=0.35;
A_1 = -29;
A_2 = -4 * A_1;

alpha_1 = 0.15;
alpha_2 = 0.5;

dydt(3) = A_1*dlambda_dt-alpha_1*Q1;
dydt(4) = A_2*dlambda_dt-alpha_2*Q2;

if (Q < 0.0)
  Qfac = (a*Q+1.0)/(1.0-Q);
else
  Qfac = (1.0+(a+2.0)*Q)/(1.0+Q);
end

%-------------------------------------------------------------------------------
% Force
%-------------------------------------------------------------------------------

Lfac = max(0, 1 + beta_0 * (lambda + min(0.87,lambda) - 1.87) );
Tension = Qfac * Tref * Lfac * XB;
