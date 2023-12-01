function varargout=ttestcorr(alpha,N,lag,delta,makefigure)
% The program is used  to determine the critical values  of
% cross-correlation through t-tests at the alpha significance level.
%
% INPUT:
% alpha             The significance level
% N                   The length of datasets 
% lag                 The time shift 
% delta              The correction factor for degrees of freedom
% makefigure    Draw  figure or not
%
% OUTPUT:
%  rc                  The critical values  of cross-correlation
%  edof              The effective dgrees of freedom
%
% Last modified by Chunchun GAO  at SDUST, 2023.09.14
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
defval('alpha',[0.1 0.05 0.01 0.001]) 
defval('N',502) 
defval('lag',(3-N):(N-3))  
defval('delta',1); 
defval('makefigure',0) 

edof=((N-abs(lag)).*delta)-2;   %Equation (4)
edof(edof<1)=1;
nalpha=length(alpha);
for nn=1:nalpha
tpv=tinv(1-alpha(nn)/2,edof);
rc(nn,:)=(N-abs(lag))/N.*tpv./sqrt(edof+tpv.^2);   %Equation (3)
end

rc=rc';

if makefigure==1
   plot(lag',rc)
   ylabel('Correlation Coefficient');xlabel('Time Shift')
   grid on
   set(gca,'GridLineStyle','--')
end

edof=edof';

 varns={rc,edof};
 varargout=varns(1:nargout);