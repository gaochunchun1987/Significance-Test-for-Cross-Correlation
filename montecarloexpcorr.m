function varargout=montecarloexpcorr(nts,confid,nexp,ntest,forcenew,datatype,data,makefigure,datafile,lagsrange)
%The program is used to determine the critical values  of
% cross-correlation by Monte Carlo experiments at the "confid" confidence level.
%
% INPUT:
% nts                    The length of time series
% confid               The confidence level
% nexp                 Monte Carlo test times
% ntest                The number of repetitions of Monte Carlo experiment
% forcenew         Whether you want to force a new generation of a save file or not 
% datatype          The signal type will be simulated by Monte Carlo
% data                 The data represents the AR1 coefficients, e.g., [0.8 0.6], if the datetype is rednoise; 
%                          or it represents the input time series if the datetype is inputdata
% makefigure      Draw  figure or not
% datafile            The file path for data storage
% lagsrange        The range of time shifts, e.g., 0:12
%
% OUTPUT:
% rtscvmean                  The critical values  of cross-correlation at time shift equal to 0 
% xrtscvmean                The critical values  of cross-correlation at certain time shift point
% xrtsmaxcvmean          The critical values  of cross-correlation within  all symmetrical time shift range
% xrtsallcvmean             The critical values  of cross-correlation for all time shift points
% xrtsmaxcvspmean      The critical values  of cross-correlation within  any specified time shift range
%
% Last modified by Chunchun GAO  at SDUST, 2023.09.14
% Last modified by Chunchun GAO  at SDUST, 2023.11.13
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
defval('nts',502) 
defval('confid',0.95) 
defval('nexp',10000) 
defval('ntest',1)  
defval('forcenew',0) 
defval('datatype','whitenoise') 
defval('data',[])
defval('makefigure',0) 
defval('datafile','d:/matlab/correlation') 
defval('lagsrange',[])

if isempty(lagsrange)
  if strcmp(datatype,'whitenoise')
    fnpl=sprintf('%s/monte_carlo_exp_crosscorr_%d_%d_%d_%f_%s.mat',datafile,nts,nexp,ntest,prod(confid),datatype);
  elseif strcmp(datatype,'rednoise')
    fnpl=sprintf('%s/monte_carlo_exp_crosscorr_%d_%d_%d_%f_%s_%f_%f.mat',datafile,nts,nexp,ntest,prod(confid),datatype,data(1),data(2));
  elseif strcmp(datatype,'inputdata')
    fnpl=sprintf('%s/monte_carlo_exp_crosscorr_%d_%d_%d_%f_%s_%f.mat',datafile,nts,nexp,ntest,prod(confid),datatype,sum(abs(data(:))));
  end
else
    if strcmp(datatype,'whitenoise')
     fnpl=sprintf('%s/monte_carlo_exp_crosscorr_%d_%d_%d_%f_%s_%d_%d.mat',datafile,nts,nexp,ntest,prod(confid),datatype,lagsrange(1),lagsrange(end));
   elseif strcmp(datatype,'rednoise')
     fnpl=sprintf('%s/monte_carlo_exp_crosscorr_%d_%d_%d_%f_%s_%f_%f_%d_%d.mat',datafile,nts,nexp,ntest,prod(confid),datatype,data(1),data(2),lagsrange(1),lagsrange(end));
   elseif strcmp(datatype,'inputdata')
     fnpl=sprintf('%s/monte_carlo_exp_crosscorr_%d_%d_%d_%f_%s_%f_%d_%d.mat',datafile,nts,nexp,ntest,prod(confid),datatype,sum(abs(data(:))),lagsrange(1),lagsrange(end));
    end
end
 

if exist(fnpl,'file')==2 && forcenew==0
    load(fnpl)
    fprintf('%s loaded by montecarloexpcorr\n',fnpl)
else
for n1=1:ntest   
    if strcmp(datatype,'whitenoise')  %Generate random white noise time series
       parfor nn=1:nexp
             ts=randn(nts,2);  
             tsnew=zscore(ts);
             xrts=crosscorr(tsnew(:,1),tsnew(:,2),'NumLags',nts-1); 
             xrtsnn(nn,:)=xrts';
       end
    elseif strcmp(datatype,'rednoise')  
        % Generate random red noise time series  by the first order autoregressive (AR1)
       parfor nn=1:nexp
              ts=[rednoise(nts,data(1)) rednoise(nts,data(2))]; 
              tsnew=zscore(ts);
              xrts=crosscorr(tsnew(:,1),tsnew(:,2),'NumLags',nts-1); 
              xrtsnn(nn,:)=xrts';
       end
    elseif strcmp(datatype,'inputdata') 
        % Generate random surrogate time series with similar autocorrelation structure to the input data
        nts=size(data,1);
        parfor nn=1:nexp
             ts=surrogatedata_fft(data,'autocorr');   
             tsnew=zscore(ts);
             xrts=crosscorr(tsnew(:,1),tsnew(:,2),'NumLags',nts-1); 
             xrtsnn(nn,:)=xrts';
        end 
    end

  %Sort and get the value of the corresponding confidence level
  xrtsdescend=sort(abs(xrtsnn));
  xrtscv(n1,:,:)=xrtsdescend(round(nexp*confid),:);
  xrtsnnabs=abs(xrtsnn)';
  for n2=1:nts-1
    xrtsnnabslag=xrtsnnabs(nts-n2:nts+n2,:);
    xrtsmaxdescend=sort(max(xrtsnnabslag));
    xrtsmaxcv(n1,n2,:)=xrtsmaxdescend(round(nexp*confid));
    xrtsuballdescend=sort(xrtsnnabslag(:))';
    xrtsallcv(n1,n2,:)=xrtsuballdescend(round(length(xrtsuballdescend)*confid));
  end
  if ~isempty(lagsrange)
      ts=randn(nts,2);  
       tsnew=zscore(ts);
       [~,lags]=crosscorr(tsnew(:,1),tsnew(:,2),'NumLags',nts-1); 
       pstart=find(lags==(lagsrange(1)));    %  Find the starting position of the time shift range in lags
       pend=find(lags==(lagsrange(end))); %  Find the final position of the time shift range in lags
       xrtsnnabslagsp=xrtsnnabs(pstart:pend,:);
      xrtsmaxdescendsp=sort(max(xrtsnnabslagsp));
      xrtsmaxcvsp(n1,:)=xrtsmaxdescendsp(round(nexp*confid));
  end    
end

% if ntest>1, average them
if ntest==1
    xrtscvmean=squeeze(xrtscv)';
else
    xrtscvmean=squeeze(mean(xrtscv))';
end

if size(xrtscvmean,1)==1
    xrtscvmean=xrtscvmean';
end
rtscvmean=xrtscvmean(nts,:);

if ntest==1
    xrtsmaxcvmean=squeeze(xrtsmaxcv);
else
   xrtsmaxcvmean=squeeze(mean(xrtsmaxcv,1));
end

if size(xrtsmaxcvmean,1)==1
    xrtsmaxcvmean=xrtsmaxcvmean';
end

if ntest==1
   xrtsallcvmean=squeeze(xrtsallcv);
else
  xrtsallcvmean=squeeze(mean(xrtsallcv,1));
end

if size(xrtsallcvmean,1)==1
    xrtsallcvmean=xrtsallcvmean';
end

 if ~isempty(lagsrange)
     if ntest==1
        xrtsmaxcvspmean=xrtsmaxcvsp;
     else
        xrtsmaxcvspmean=mean(xrtsmaxcvsp);
     end
 end
 
if isempty(lagsrange)
   save(fnpl,'rtscvmean','xrtscvmean','xrtsmaxcvmean','xrtsallcvmean');
else
   save(fnpl,'rtscvmean','xrtscvmean','xrtsmaxcvmean','xrtsallcvmean','xrtsmaxcvspmean'); 
end

end

if makefigure==1
     if strcmp(datatype,'whitenoise')
          ts=randn(nts,2);  
    elseif strcmp(datatype,'rednoise')
          ts=[rednoise(nts,data(1)) rednoise(nts,data(2))]; 
    elseif strcmp(datatype,'inputdata')
         ts=data; 
    end
      tsnew=zscore(ts);
      xrts=crosscorr(tsnew(:,1),tsnew(:,2),'NumLags',nts-1); 
      lags=(1-nts:nts-1)';
      plot(lags,xrts,'b');title('Cross Correlation');ylabel('Coefficient');xlabel('Lag')
      hold on
      plot(lags,[repmat(rtscvmean,size(lags)) repmat(rtscvmean*-1,size(lags))],'--k'...
        ,lags,[repmat(xrtsmaxcvmean(end,:),size(lags)) repmat(xrtsmaxcvmean(end,:)*-1,size(lags))],'--r')
      plot(lags,[repmat(xrtsallcvmean(end,:),size(lags)) repmat(xrtsallcvmean(end,:)*-1,size(lags))],'--g'...
        ,lags,[xrtscvmean xrtscvmean*-1],'--m')
      text([350 350],[rtscvmean(1)+0.01 rtscvmean(1)*-1-0.01],'95%(Lag=0)');
%   text([350 350],[rtscvmean(2)+0.01 rtscvmean(2)*-1-0.01],'99%');
      text([350 350],[xrtsmaxcvmean(end,1)+0.01 xrtsmaxcvmean(end,1)*-1-0.01],'95%(Max)','Color','r');
%   text([350 350],[xrtsmaxcvmean(end,2)+0.01 xrtsmaxcvmean(end,2)*-1-0.01],'99%(Max)','Color','r');
      text([350 350],[xrtsallcvmean(end,1)+0.01 xrtsallcvmean(end,1)*-1-0.01],'95%(All)','Color','g');
      text([350 350],[xrtscvmean(end,1)+0.05 xrtscvmean(end,1)*-1-0.05],'95%(Each)','Color','m');
end

if isempty(lagsrange)
    varns={rtscvmean,xrtscvmean,xrtsmaxcvmean,xrtsallcvmean};
else
    varns={rtscvmean,xrtscvmean,xrtsmaxcvmean,xrtsallcvmean,xrtsmaxcvspmean};
end
varargout=varns(1:nargout);