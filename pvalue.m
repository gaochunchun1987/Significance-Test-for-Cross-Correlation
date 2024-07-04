function varargout=pvalue(XX,YY,r,lag,lagsrange,edofm)
% The program is used  to determine the P-values  of
% cross-correlation through t-tests.
%
% INPUT:
% XX                  The time series X
% YY                  The time series Y
% r                     The correlation coefficients
% lag                 The time shift 
% lagsrange      The time shifts range
% delta              The correction factor for degrees of freedom
%
% OUTPUT:
%  pv                  The P-values
%
% Last modified by Chunchun GAO  at SDUST, 2024.06.26
% Email: gaochunchun@sdust.edu.cn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CAUTION: THE SOFTWARE AND ITS ALGORITHMS ARE EXCLUSIVELY AVAILABLE FOR INDIVIDUAL 
% USERS TO ACQUIRE KNOWLEDGE AND EMPLOY IN SCIENTIFIC  RESEARCH. IT IS STRICTLY
% PROHIBITED FOR ANY USER TO EXPLOIT THE SOFTWARE AND ALGORITHMS FOR COMMERCIAL
% PURPOSES (INCLUDING, BUT NOT LIMITED TO,  EMPLOYING THE SOFTWARE IN GOVERNMENT
% PROCUREMENT OR BIDDING PROCESSES).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------------------------
% Set the default values of the input variables
defval('XX',rednoise(502,0.8))
defval('YY',rednoise(502,0.8))
N=length(XX);
defval('r',0.325)
defval('lag',0)  
defval('lagsrange',0)
defval('edofm','xBH')

delta=edofcf(XX,YY,lag,edofm,0); % Equations 4
edof=((N-abs(lag)).*delta)-2;   %Equation (4)
edof(edof<1)=1;

rv=r/((N-abs(lag))/N);
v=edof;

%TPVALUE Compute p-value for t statistic.
Tstat = rv .* sqrt(v./ (1 - rv.^2));
x=-abs(Tstat);

normcutoff = 1e7;
if length(x)~=1 && length(v)==1
   v = repmat(v,size(x));
end

% Initialize P.
p = NaN(size(x), 'like', x([])); % MOD 1
nans = (isnan(x) | ~(0<v)); % v == NaN ==> (0<v) == false

% First compute F(-|x|).
%
% Cauchy distribution.  See Devroye pages 29 and 450.
cauchy = (v == 1);
p(cauchy) = .5 + atan(x(cauchy))/pi;

% Normal Approximation.
normal = (v > normcutoff);
p(normal) = 0.5 * erfc(-x(normal) ./ sqrt(2));

% See Abramowitz and Stegun, formulas 26.5.27 and 26.7.1.
gen = ~(cauchy | normal | nans);
p(gen) = betainc(v(gen) ./ (v(gen) + x(gen).^2), v(gen)/2, 0.5)/2;

% Adjust for x>0.  Right now p<0.5, so this is numerically safe.
reflect = gen & (x > 0);
p(reflect) = 1 - p(reflect);

% Make the result exact for the median.
p(x == 0 & ~nans) = 0.5;
p=2*p;

if strcmp(edofm,'xBH')
   edofm='BH';
end
   deltatau=edofcf(XX,YY,lagsrange,edofm,0);
  %The number of independent time shift points within the range   Equation(6)
   ni=deltatau*length(lagsrange);  
if ni<1
     ni=1;
end
 
pv=1-(1-p).^ni; 

 varns={pv};
 varargout=varns(1:nargout);