%% The program requires access to the path "d:/matlab/correlation", so the corresponding file name must be established on disk 'd' before running
%% Last modified by Chunchun GAO  at SDUST, 2023.11.15
clear;
close all
N=500;

% Test for white noise
% Two random white noise signals are generated
XX=zscore(randn(N,1));
YY=zscore(randn(N,1));
confid=0.95;
alpha=1-confid;
lagsrange=(N*0.5)*-1:(N*0.5);
% The significance testing for cross correlation is conducted at the known time shift point
[rc,edof,signresult,signxrts]=cvttest(XX,YY,lagsrange,'WN',alpha,1,0); % Theoretical method 
[~,rcmc]=montecarloexpcorr(N,confid,100000,1,0,'whitenoise',[],0,[],[]); % Monte Carlo method
%  The significance testing for cross correlation within a specified range of time shifts
[rcrg,edofrg,signresultrg,signxrtsrg]=cvttest(XX,YY,lagsrange,'WN',alpha,length(lagsrange),0);
[~,~,~,~,rcrgmc]=montecarloexpcorr(N,confid,100000,1,0,'whitenoise',[],0,[],lagsrange); % Monte Carlo method
 
 % Test for red noise
 % Two random red noise signals are generated
XXred=zscore(rednoise(N,0.8));
YYred=zscore(rednoise(N,0.8));
% The significance testing for cross correlation is conducted at the known time shift point
[rcred,edofred,signresultred,signxrtsred]=cvttest(XXred,YYred,lagsrange,'xBH',alpha,1,0); % Theoretical method 
[~,rcmcred]=montecarloexpcorr(N,confid,100000,1,0,'inputdata',[XXred,YYred],0,[],[]); % Monte Carlo method
%  The significance testing for cross correlation within a specified range of time shifts
deltatau=edofcf(XXred,YYred,lagsrange,'BH',0);
ni=deltatau*length(lagsrange);  %The number of independent time shift points within the range   Equation(6)
[rcrgred,edofrgred,signresultrgred,signxrtsrgred]=cvttest(XXred,YYred,lagsrange,'xBH',alpha,ni,0);
[~,~,~,~,rcrgmcred]=montecarloexpcorr(N,confid,100000,1,0,'inputdata',[XXred,YYred],0,[],lagsrange); % Monte Carlo method

% plotting
figure(1)
clf
subplot(2,2,1)
[xrts,lagscc]=crosscorr(XX,YY,'NumLags',N-1);  
plot(lagscc,xrts)
hold on 
 plot(lagsrange',[rc -1*rc],'go')
 pstart=find(lagscc==(lagsrange(1)));    %  Find the starting position of the time shift range in lagscc
 pend=find(lagscc==(lagsrange(end))); %  Find the end position of the time shift range in lagscc
 plot(lagsrange',[rcmc(pstart:pend,:) -1*rcmc(pstart:pend,:)],'r-')
 ylabel('Correlation Coefficient');xlabel('Time Shift')
 grid on
 set(gca,'GridLineStyle','--')
 title('(a) Time Shift Point (White Noise)')
subplot(2,2,2)
plot(lagscc,xrts)
  hold on 
 plot(lagsrange',[rcrg -1*rcrg],'go')
 plot(lagsrange',[repmat(rcrgmc,length(lagsrange),1) -1*repmat(rcrgmc,length(lagsrange),1)],'r-') 
  ylabel('Correlation Coefficient');xlabel('Time Shift')
 grid on
 set(gca,'GridLineStyle','--')
 title('(b) Time Shift Range (White Noise)')
subplot(2,2,3)
[xrtsred,lagscc]=crosscorr(XXred,YYred,'NumLags',N-1);  
plot(lagscc,xrtsred)
  hold on 
 plot(lagsrange',[rcred -1*rcred],'go')
 pstart=find(lagscc==(lagsrange(1)));    %  Find the starting position of the time shift range in lagscc
 pend=find(lagscc==(lagsrange(end))); %  Find the final position of the time shift range in lagscc
 plot(lagsrange',[rcmcred(pstart:pend,:) -1*rcmcred(pstart:pend,:)],'r-')
 ylabel('Correlation Coefficient');xlabel('Time Shift')
 grid on
 set(gca,'GridLineStyle','--')
  title('(c) Time Shift Point (Red Noise)')
subplot(2,2,4)
plot(lagscc,xrtsred)
  hold on 
 plot(lagsrange',[rcrgred -1*rcrgred],'go')
 plot(lagsrange',[repmat(rcrgmcred,length(lagsrange),1) -1*repmat(rcrgmcred,length(lagsrange),1)],'r-') 
  ylabel('Correlation Coefficient');xlabel('Time Shift')
 grid on
 set(gca,'GridLineStyle','--')
  title('(d) Time Shift Range (Red Noise)')

