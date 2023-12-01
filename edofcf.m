function varargout=edofcf(XX,YY,lag,edofm,makefigure)
% The program is utilized for computing the correction factors of effective degrees of freedom, 
% with these methodologies being summarized in Afyouni et al. (2019)
%
% INPUT:
% XX                 The time series X
% YY                  The time series Y
% lag                 The time shifts tau
% edofm            The method for computing correction factors for degrees of freedom
% makefigure    Draw  figure or not
%
% OUTPUT:
% delta              The correction factor for degrees of freedom
%
% Last modified by Chunchun GAO  at SDUST, 2023.09.20
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
defval('XX',rednoise(100,0.8))
defval('YY',rednoise(100,0.6))
N=length(XX);
defval('lag',5-N:N-5)
defval('edofm','xBH')
defval('makefigure',1)
acx = autocorr(XX,N-1);
acy = autocorr(YY,N-1);
 
if strcmp(edofm,'WN') % white noise
    delta=1;
elseif strcmp(edofm,'B35')  % Bartlett,1935
    acxy1 =acx(2)*acy(2);
    delta= (1-acxy1)./(1+acxy1);
elseif strcmp(edofm,'Q47')  % Quenouille, 1947
    acxy =acx(2:end).*acy(2:end);
    delta= 1/(1+2*sum(acxy));
elseif strcmp(edofm,'GQ47')  % Fox et al., 2005
    acxy =(acx(2:end).^2+acy(2:end).^2)/2;
    delta= 1/(1+2*sum(acxy));   
elseif strcmp(edofm,'BH')  % Bayley and Hammersley,1946
    wgt    = (N-1:-1:1)'./N;
    acxy =wgt.*acx(2:end).*acy(2:end);
    delta= 1/(1+2*sum(acxy));
elseif strcmp(edofm,'xDF')  %Afyouni et al., 2019
    v_xdf=xDF([XX YY]',N,'truncate','adaptive','TVOn');
    r=corr(XX, YY);
    delta=  (1-r^2)^2/ v_xdf(1,2)/N;   
elseif strcmp(edofm,'xBH')  % Equation (5) in Gao et al.,2023  
    for nn=1:length(lag)
       wgt = (N-abs(lag(nn))-1:-1:1)'./(N-abs(lag(nn)));
       acxy =wgt.*acx(2:end-abs(lag(nn))).*acy(2:end-abs(lag(nn)));
       delta(nn)= 1/(1+2*sum(acxy));
    end   
end

if any(any(delta>1)) 
    delta(delta>1) = 1;
end

if makefigure==1
plot(lag,delta)
end

 varns={delta};
 varargout=varns(1:nargout);