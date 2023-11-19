 function varargout=cvttest(XX,YY,lagsrange,edofm,alpha,ni,makefigure)
% The program is utilized for determining the critical values of cross-correlation through t-tests 
% at the adjusted significance level, ensuring statistical significance in hypothesis testing.
%
% INPUT:
% XX                  The time series X
% YY                  The time series Y
% lagsrange      The time shifts range
% edofm            The method for computing correction factors for degrees of freedom, see "edofcf" in detail
% alpha              The significance level      
% ni                   The number of independent time shift points within the  range
% makefigure    Draw  figure or not
%
% OUTPUT:
%  rc                   The critical values  of cross-correlation
%  edof               The effective dgrees of freedom
% signresult       Significant (1) or not (0)
% signxrts          Significant cross correlation and time shifts
%
% Last modified by Chunchun GAO  at SDUST, 2023.09.28
%-------------------------------------------------------------------------------------------
% Set the default values of the input variables
defval('XX',randn(500,1))
defval('YY',randn(500,1))
N=length(XX);
defval('lagsrange',round(N*0.5)*-1:round(N*0.5))
defval('edofm','WN')
defval('alpha',0.05) 
defval('ni',length(lagsrange)); 
if size(lagsrange,2)==1
    lagsrange=lagsrange';
end
defval('makefigure',1)

 delta=edofcf(XX,YY,lagsrange,edofm,0); % Equations 4
 % compute the cross correlation
 [xrts,lagscc]=crosscorr(XX,YY,'NumLags',N-1);  
pstart=find(lagscc==(lagsrange(1)));    %  Find the starting position of the time shift range in lagscc
pend=find(lagscc==(lagsrange(end))); %  Find the final position of the time shift range in lagscc
xrtscmp=xrts(pstart:pend,1);
 
 if ni<1
     ni=1;
 end
 
 alphac=1-(1-alpha).^(1/ni); % Equation 7
 [rc,edof]=ttestcorr(alphac,N,lagsrange,delta); % Equation 8

for nn=1:length(alpha)
signtest=find(abs(xrtscmp)>rc(:,nn)); % Determine whether it is significant
if isempty(signtest)
    signresult(nn)=0;
else
    signresult(nn)=1;
end
signxrts{nn}=[xrtscmp(signtest) lagsrange(signtest)'];
end
 
 if makefigure==1
     plot(lagscc,xrts)
     hold on 
     plot(lagsrange',[rc -1*rc],'r.')
     hold off
     ylabel('Correlation Coefficient');xlabel('Time Shift')
     grid on
     set(gca,'GridLineStyle','--')
 end
 
 varns={rc,edof,signresult,signxrts};
 varargout=varns(1:nargout);