function varargout=surrogatedata_fft(data,gmethod,makefigure)
% The program is used to generate random surrogate data according to Prichard and Theiler(1994)
%
% INPUT:
% data                   Input data
% gmethod           The method of  generate random surrogate data
% makefigure        Draw figure or not
%
% OUTPUT:
% sdata                  The  surrogate data
%
% Last modified by Chunchun GAO  at SDUST, 2023.06.17
%-------------------------------------------------------------------------------------------
% Set the default values of the input variables
defval('data',[rednoise(500,0.8) rednoise(500,0.6)])
defval('gmethod','autocorr')
defval('makefigure',0)

dataf=fft(data);
dataang=angle(dataf);
  %Generate data with similar autocorrelation structure
if strcmp(gmethod,'autocorr')
randang=rand(size(dataf))*2*pi;
%Generate data with similar autocorrelation and cross-correlation structures
elseif strcmp(gmethod,'crosscorr') 
 randang1=rand(size(dataf(:,1)))*2*pi;
 randang=repmat(randang1,[1 size(dataf,2)]);
end
sdataf=abs(dataf).*exp(i*(angle(dataf)+randang));   
sdata=real(ifft(sdataf));

if makefigure==1
nts=size(data,1);
x=data(:,1);
y=data(:,2);
x1=sdata(:,1);
y1=sdata(:,2);
xrts=crosscorr(x,y,'NumLags',nts-1);
xrts1=crosscorr(x1,y1,'NumLags',nts-1);
 lags=(1-nts:nts-1)'; 
 acx=AC_fft(x,nts);
 acy=AC_fft(y,nts);
 acx1=AC_fft(x1,nts);
 acy1=AC_fft(y1,nts);
subplot(2,3,1),plot((1:nts)',[x x1])
subplot(2,3,2),plot((1:nts)',[y y1])
subplot(2,3,3),plot(0:nts-1,[acx; acx1])
subplot(2,3,4),plot(0:nts-1,[acy; acy1])
subplot(2,3,5),plot(lags,[xrts xrts1])
end

varns={sdata};
varargout=varns(1:nargout);