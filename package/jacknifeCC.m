function [fitI,info1,I,Q]=jacknifeCC(X,K,samples,varargin)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% JACKNIFECC extrapolates channel capacity for a set of (multidimensional) inputs,
% X. Actual information calculation is performed by "getCC".
%
%
% X         Cell array (size=nx1) of sets of individual responses to n different inputs
% K         Kth nearest neighbor (sole parameter passed to KNN algorithm)
% samples   Number of subsets to measure before computing jacknife extrapolation. Subsets 
%           are randomly drawn to contain between 65%-90% of total samples. The number of 
%           "draws" per sample is computed as 5/(x%)^2.
% varargin  If initialized, forces "verbose mode" (text+graphical output)
%
%
% fitI      Extrapolated mutual information for entire set
% info1     Mutual information for full set
% I         Mutual information for all subsets
% Q         "Optimal" weights (of all possible inputs) that maximizes mutual information
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if nargin<4
    verbose = 0;
else
    verbose = 1;
end

Nsig=length(X);

[info1, Q] = getCC_Previous(X,K,0); %full set of data 

if verbose
 % #bins, #samples, #draws for jackknife to correct for sample bias        
SampleSize=linspace(0.95,0.6,samples); %fraction of total        
D= round(5./SampleSize(1:end).^2); %# of draws for each # of samples

Nss=samples; %# of subsets for jackknife

Ik=cell(Nss,1);
for k=1:Nss 
    if verbose; fprintf(['subset' num2str(k)]); end
    info=zeros(D(k),1);
    for ds=1:D(k)    
        if verbose; fprintf('.'); end
        sX=cell(Nsig,1);
        for i=1:Nsig               
            rp=randperm(size(X{i},2));            
            sX{i}=X{i}(:,rp(1:round(SampleSize(k)*size(X{i},2))));                  
        end
        info(ds) = getCC_Previous(sX,K,0);   
    end
    Ik{k}=info;
end
I=[info1; Ik];
if verbose; fprintf('\n'); end

meanI=zeros(size(I));
stdI=zeros(size(I));
for i=1:length(I)
    meanI(i)=mean(I{i});
    stdI(i)=std(I{i});
end
ss=[1 SampleSize]'; 
pI=polyfit(1./ss,meanI,1);
fitI=pI(2);

if verbose
   fig1=figure('position',[50 50 500 500]);
   hold on
   plot(1./ss,meanI,'.')
   plot(1,info1,'.')
   plot([0,2],pI(2)+pI(1)*[0,2],'r')
   xlabel('1/SS')
   ylabel('info (bits)')
   
   %save jacknife_fig_data ss meanI pI stdI
   
   %exportfig(gcf,'jacknife','Format','eps','Color','rgb',...
   % 'FontSizeMin',8,'FontSizeMax',12,'FontSize',2,...
   % 'LineWidthMax',1,'LineWidthMin',0.25)
end

else
    fitI=0;I=0;
end
