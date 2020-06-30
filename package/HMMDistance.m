%Quantification on the training goodness of hidden Markov model, including
%False nearest neighbor probability, Relative KL-divergence, Rescaled likelihood.

function HMM=HMMDistance(estTR,estE,OtherPara,seqs,Y,Nq,foldername,ID,SampledTraj)


DistanceConditions=0;
    TrajSample=cell(1,Nq);
    TrajSample_Smooth=cell(1,Nq);
    Y_Smooth=cell(1,Nq);
    seqs_Smooth=cell(1,Nq);
    TrajSample_HiddenStates=cell(1,Nq);
    estTRsort=cell(1,Nq);
    estEsort=cell(1,Nq);
    DistanceTrajs=cell(1,Nq);Distance2_Trajs=cell(1,Nq);
    
    
   nbins=10;


IniDist=cell(Nq,1);
count1=0;count2=0;
for Id_index=1:Nq%+1
     seqs{Id_index}=seqs{Id_index}*OtherPara.Conversion;
disp(Id_index);
if OtherPara.Dataset==1
if min(abs(ID(Id_index)-[548,546,778,780,720,610,779]))==0
        YUpLimit=5;
    elseif min(abs(ID(Id_index)-[650,566]))==0
        YUpLimit=4;
    elseif min(abs(ID(Id_index)-[756]))==0
        YUpLimit=7;
    else
        YUpLimit=OtherPara.MaxValue;
end
    YUpLimit2=OtherPara.state;
    YUpLimit=5;%max(max(Y{Id_index}));% 7 is suggested by their paper...
    YUpLimit=OtherPara.YUpLimit(Id_index);
else
    YUpLimit2=OtherPara.state;
    %if Id_index==1
    YUpLimit=OtherPara.YUpLimit(Id_index);%OtherPara.MaxValue;
    %elseif Id_index==2
     %   YUpLimit=OtherPara.MaxValue;
    %end
end
  
    %assignin('base', 'YUpLimit', YUpLimit);

    sampleSize=round(OtherPara.NumTraj(Id_index)*OtherPara.EnhancedSample);
    if sampleSize<5
        sampleSize=500;
    end
    if OtherPara.NumTraj(Id_index)<5
        OtherPara.NumTraj(Id_index)=500;
    end
    TrajLength=OtherPara.TrajLength(Id_index); %*2 is to predict the sampling for longer time

opt= maketicks(1:TrajLength,[0, YUpLimit],0); opt.Name = 'ddd';
opt2= maketicks(1:TrajLength,[0, YUpLimit2],0); opt2.Name = 'ddd';


HM= figure ('position', [00, 10, 600, 800]);%('name','650');  
    %assignin('base', 'TrajSample_Smooth', TrajSample_Smooth);
     for jjj=1:size(seqs{Id_index},1)
    seqs_Smooth{1,Id_index}(jjj,:) = max(smooth(1:size(seqs{Id_index},2),...
    seqs{Id_index}(jjj,:),0.07,'loess'),0);
    end
     [metrics,aux] = nfkbmetrics(seqs_Smooth{Id_index});%(:,1:OtherPara.TrajLength(Id_index)));
[aaa,index]=sortrows(metrics.peakfreq);
tblB=seqs_Smooth{Id_index}(flipud(index),:);  
h=colormapStack(tblB,[],opt, HM);
 hh = colorbar;  
    s.EdgeColor = 'none';view(2);
    xlabel('Time (h)');ylabel('Single cells');
    if OtherPara.Dataset==2 % 2 is p53 after Mdmx;
         xticks([1:20:122]);xticklabels({'0','10','20','30','40','50','60'});
         %xticks([1:24:130]);xticklabels({'0','2','4','6','8','10','12'});
    elseif OtherPara.Dataset==3 % 3 is Erk
         xticks([1:40:130]);xticklabels({'0','5','10','15'});
         %xticks([1:24:130]);xticklabels({'0','2','4','6','8','10','12'});
    else%elseif OtherPara.Dataset==1
        xticks([1:24:130]);xticklabels({'0','2','4','6','8','10','12'});
   
        xlim([1 max(OtherPara.TimePointsToUse)])
        TimeUnit=round(OtherPara.TotalTimeLength/6);
        xTime=0:TimeUnit:5*TimeUnit;
        xTimePoints=round(xTime*OtherPara.TimePointsToUse/OtherPara.TotalTimeLength);    
        xTimePoints(1)=1;
        xticks(xTimePoints);xticklabels(cellstr(string(xTime)));%xticklabels({'0','2','4','6','8','10','12'});
    end
   
    %xticks([1:24:130]);xticklabels({'0','2','4','6','8','10','12'});
    set(gca,'FontSize',36);
     figurenamehmm=[foldername,'\BinedData_',num2str(ID(Id_index)),'States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize),'.jpg'];
 saveas(gcf,figurenamehmm); 
 
 
 
 
 HM= figure ('position', [00, 10, 600, 800]);%('name','650');  
    %assignin('base', 'TrajSample_Smooth', TrajSample_Smooth);
     for jjj=1:size(Y{Id_index},2)
    Y_Smooth{1,Id_index}(jjj,:) = max(smooth(1:size(Y{Id_index},1),...
    Y{Id_index}(:,jjj)',0.07,'loess'),0);
    end
     [metrics,aux] = nfkbmetrics(Y_Smooth{1,Id_index});%(:,1:OtherPara.TrajLength(Id_index)));   
[aaa,index]=sortrows(metrics.peakfreq);
tblB=Y_Smooth{1,Id_index}(flipud(index),:);  
h=colormapStack(tblB,[],opt, HM);
 hh = colorbar;  
    s.EdgeColor = 'none';view(2);
    xlabel('Time (h)');ylabel('Single cells');
    if OtherPara.Dataset==2 % 2 is p53 after Mdmx;
         xticks([1:20:122]);xticklabels({'0','10','20','30','40','50','60'});
         %xticks([1:24:130]);xticklabels({'0','2','4','6','8','10','12'});
    elseif OtherPara.Dataset==3 % 3 is Erk
         xticks([1:40:130]);xticklabels({'0','5','10','15'});
         %xticks([1:24:130]);xticklabels({'0','2','4','6','8','10','12'});
    else% elseif OtherPara.Dataset==1
        xticks([1:24:130]);xticklabels({'0','2','4','6','8','10','12'});
   
        xlim([1 max(OtherPara.TimePointsToUse)])
        TimeUnit=round(OtherPara.TotalTimeLength/6);
        xTime=0:TimeUnit:5*TimeUnit;
        xTimePoints=round(xTime*OtherPara.TimePointsToUse/OtherPara.TotalTimeLength);    
        xTimePoints(1)=1;
        xticks(xTimePoints);xticklabels(cellstr(string(xTime)));%xticklabels({'0','2','4','6','8','10','12'});
    end
%    xticks([1:24:130]);xticklabels({'0','2','4','6','8','10','12'});
    set(gca,'FontSize',36);
     figurenamehmm=[foldername,'\Data_',num2str(ID(Id_index)),'States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize),'.jpg'];
 saveas(gcf,figurenamehmm); 
 
 
 
   %Sample more times to get closer trajectories
    distince=zeros(1,sampleSize);TrajSample_Temp=[];Index=[];TrajSample_HiddenStatesTemp=[];%(jjj,:)
    for jjj=1:sampleSize
       % disp(jjj);
    [seqSample,states] = hmmgenerate(TrajLength,estTR{1,Id_index},estE{1,Id_index});
 
    TrajSample_Temp(jjj,:)=seqSample*OtherPara.Conversion;%[TrajSample{1,Id_index}',seq'*OtherPara.Conversion]';%*OtherPara.Conversion
    TrajSample_HiddenStatesTemp(jjj,:)=states;
    % TrajSample_SmoothTemp(jjj,:) = max(smooth(1:size(TrajSample{1,Id_index},2),...
   % TrajSample{1,Id_index}(jjj,:),0.07,'loess'),0);
    
    %disp(size(Y{Id_index})); disp(size(TrajSample_Smooth{1,Id_index}(jjj,:)));
    [index22(jjj),DistK]=annquery(Y{Id_index},TrajSample_Temp(jjj,1:OtherPara.TrajLength(Id_index))',1); %the first is data as reference
    distince(jjj)=DistK;
    end
   % %assignin('base', 'index22', sort(index22));ddd
    
    [B,Index] = sort(distince);
    
    TrajSample_Smooth{1,Id_index}(1:OtherPara.NumTraj(Id_index),:)=TrajSample_Temp(Index(1:OtherPara.NumTraj(Id_index)),:);
    TrajSample_HiddenStates{1,Id_index}(1:OtherPara.NumTraj(Id_index),:)=TrajSample_HiddenStatesTemp(Index(1:OtherPara.NumTraj(Id_index)),:);
    
    for jjj=1:OtherPara.NumTraj(Id_index)
    TrajSample_Smooth{1,Id_index}(jjj,:) = max(smooth(1:size(TrajSample_Smooth{1,Id_index},2),...
    TrajSample_Smooth{1,Id_index}(jjj,:),0.07,'loess'),0);
    TrajSample_HiddenStates{1,Id_index}(jjj,:) = max(smooth(1:size(TrajSample_HiddenStates{1,Id_index},2),...
    TrajSample_HiddenStates{1,Id_index}(jjj,:),0.07,'loess'),0);
    end
  
    HM= figure ('position', [00, 10, 600, 800]);%('name','650');  
    %assignin('base', 'TrajSample_Smooth', TrajSample_Smooth);
     [metrics,aux] = nfkbmetrics(TrajSample_Smooth{1,Id_index});%(:,1:OtherPara.TrajLength(Id_index)));
[aaa,index]=sortrows(metrics.peakfreq);
tblB=TrajSample_Smooth{1,Id_index}(flipud(index),:);  
h=colormapStack(tblB,[],opt, HM);
 hh = colorbar;  
    s.EdgeColor = 'none';view(2);
    xlabel('Time (h)');ylabel('Single cells');
    if OtherPara.Dataset==2 % 2 is p53 after Mdmx;
         xticks([1:20:122]);xticklabels({'0','10','20','30','40','50','60'});
         %xticks([1:24:130]);xticklabels({'0','2','4','6','8','10','12'});
    elseif OtherPara.Dataset==3 % 3 is Erk
         xticks([1:40:130]);xticklabels({'0','5','10','15'});
         %xticks([1:24:130]);xticklabels({'0','2','4','6','8','10','12'});
    else%  elseif OtherPara.Dataset==1
        xticks([1:24:130]);xticklabels({'0','2','4','6','8','10','12'});
    
        xlim([1 max(OtherPara.TimePointsToUse)])
        TimeUnit=round(OtherPara.TotalTimeLength/6);
        xTime=0:TimeUnit:5*TimeUnit;
        xTimePoints=round(xTime*OtherPara.TimePointsToUse/OtherPara.TotalTimeLength);    
        xTimePoints(1)=1;
        xticks(xTimePoints);xticklabels(cellstr(string(xTime)));%xticklabels({'0','2','4','6','8','10','12'});
    end
    set(gca,'FontSize',36);
     figurenamehmm=[foldername,'\HMM_Sampling_X',num2str(OtherPara.EnhancedSample),'_',num2str(ID(Id_index)),'States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize),'.jpg'];
 saveas(gcf,figurenamehmm); 
    
    
 %%plot hidden states heat map
%    HM2= figure ('position', [00, 10, 600, 800]);%('name','650');    
% tblB=TrajSample_HiddenStates{1,Id_index}(flipud(index),:);  
% h=colormapStack(tblB,[],opt2, HM2);
%  hh = colorbar;  
%     s.EdgeColor = 'none';view(2);
%     xlabel('Time (h)');ylabel('Single cells');
%     xticks([1:24:130]);xticklabels({'0','2','4','6','8','10','12'});
%     set(gca,'FontSize',36);
%      figurenamehmm=[foldername,'\HMM_HiddenStates_FitToAll_X',num2str(OtherPara.EnhancedSample),'_',num2str(ID(Id_index)),'States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize),'.jpg'];
%  %%%%saveas(gcf,figurenamehmm); 
 
%%
%Get confusion of kNN between sampled set and data...
WholeTraj=[Y_Smooth{1,Id_index};TrajSample_Smooth{1,Id_index}];
HMM.kNNnum=60;

NumofTraj=size(TrajSample_Smooth{1,Id_index},1);kNNSampletoData=[];kNNDatatoSample=[];
%data to the whole sets
for ll=1:NumofTraj
[IndexK,DistK]=annquery(WholeTraj',TrajSample_Smooth{1,Id_index}(ll,:)',HMM.kNNnum); %the first is data as reference
kNNSampletoData(ll)= length(find(IndexK < NumofTraj));
end

%sampled set to the whole sets
% disp(OtherPara.NumTraj(Id_index))
% disp(size(Y_Smooth{1,Id_index}));
% disp(NumofTraj);
for ll=1:NumofTraj
[IndexK,DistK]=annquery(WholeTraj',Y_Smooth{1,Id_index}(ll,:)',HMM.kNNnum); %the first is data as reference
kNNDatatoSample(ll)= length(find(IndexK > NumofTraj));
end

%assignin('base', 'kNNSampletoData', kNNSampletoData);
%assignin('base', 'kNNDatatoSample', kNNDatatoSample);


 %assignin('base', 'kNNSampletoData', kNNSampletoData);
% %assignin('base', 'kNNDatatoSample', kNNDatatoSample);
%disp(NumofTraj);
 HM= figure ('position', [00, 10, 700, 1000]);
 subplot(2,1,1)
 hhh=histogram(kNNSampletoData,round(HMM.kNNnum/3));hold on;
 pd = fitdist(kNNSampletoData','Binomial','NTrials',HMM.kNNnum);
 HMM.FalseNeighborProb1(Id_index)=pd.p;
 
 x_values = [0:0.1:HMM.kNNnum];%[hhh.BinEdges(1):0.01:hhh.BinEdges(end)];
 y = binopdf(x_values,pd.N,pd.p);%pdf(pd,x_values);
 plot(x_values,y/max(y)*max(hhh.Values),'linewidth',2);
 %assignin('base', 'aap', pd);
 %assignin('base', 'ab1', x_values);
 %assignin('base', 'ab2', y);
 %assignin('base', 'aa1', hhh);
 
   %h=legend(num2str(round(HMM.FalseNeighborProb1(Id_index),2)),'Location','best');
   h=text(0.6,0.9,['Mean Prob. ',num2str(round(HMM.FalseNeighborProb1(Id_index),2))],'Units','normalized');
  set(h,'FontSize',20);
  xlabel(['Number of ',num2str(HMM.kNNnum),'-nearest neighbor from sample']);
  ylabel('Count of trajectories');
 xlim([0 HMM.kNNnum])
 title([num2str(HMM.kNNnum),'-nearest neighbor of sample to data']);
 set(gca,'FontSize',20);
  subplot(2,1,2)
 hhh=histogram(kNNDatatoSample,round(HMM.kNNnum/3));hold on;
 pd2 = fitdist(kNNDatatoSample','Binomial','NTrials',HMM.kNNnum);
 HMM.FalseNeighborProb12(Id_index)=pd2.p;
 
   x_values = [0:0.1:HMM.kNNnum];%  x_values = hhh.BinEdges;%[hhh.BinEdges(1):0.01:hhh.BinEdges(end)];
 y = binopdf(x_values,pd2.N,pd2.p);%pdf(pd,x_values);;
 plot(x_values,y/max(y)*max(hhh.Values),'linewidth',2);
 
 
  %h=legend(num2str(round(HMM.FalseNeighborProb12(Id_index),2)),'Location','best');
   h=text(0.6,0.9,['Mean Prob. ',num2str(round(HMM.FalseNeighborProb12(Id_index),2))],'Units','normalized');
  set(h,'FontSize',20);
   xlabel(['Number of ',num2str(HMM.kNNnum),'-nearest neighbor from sample']);
  ylabel('Count of trajectories');
  xlim([0 HMM.kNNnum])
 title([num2str(HMM.kNNnum),'-nearest neighbor of sample to data']);
set(gca,'FontSize',20); 
figurenamehmm=[foldername,'\kNNrelation_k',num2str(HMM.kNNnum),'_X',num2str(OtherPara.EnhancedSample),'_',num2str(ID(Id_index)),'States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize),'.jpg'];   
 %%saveas(gcf,figurenamehmm);

%  HM= figure ('position', [00, 10, 600, 800]);
%  subplot(2,1,1)
%  histogram(kNNSampletoData)
%  xlim([0 HMM.kNNnum])
%  title('kNN of sample to data');
%  set(gca,'FontSize',20);
%   subplot(2,1,2)
%  histogram(kNNDatatoSample)
%   xlim([0 HMM.kNNnum])
%  title('kNN of data to sample');
% set(gca,'FontSize',20); 
% figurenamehmm=[foldername,'\kNNrelation_k',num2str(HMM.kNNnum),'_X',num2str(OtherPara.EnhancedSample),'_',num2str(ID(Id_index)),'States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize),'.jpg'];   
%  %%saveas(gcf,figurenamehmm);


%%
 if OtherPara.DistanceAsDistribution==1
     KLDivergence=[];KLEntropy=[];KLDivergenceToData=[];KLEntropyToData=[];
     for kk=1:OtherPara.TrajLength(Id_index)
%          h1=histogram(seqs{Id_index}(:,kk),nbins);
%          h2=histogram(TrajSample_Smooth{1,Id_index}(:,kk),nbins);
%          KLDivergence(kk)=KLDiv(h1.Values,h2.Values);
          [h0,edges] = histcounts(Y{Id_index,1}(kk,:)',nbins); %for KL divergence to data
          [h1,edges] = histcounts(seqs{Id_index}(:,kk),nbins); %for KL divergence to binned seq
          [h2,edges] = histcounts(TrajSample_Smooth{1,Id_index}(:,kk),nbins);
          %assignin('base', 'h2', h2);assignin('base', 'h1', h1);
%          [KLDivergence(kk) KLEntropy(kk)]=KLDiv(h0,h2);
         [KLDivergence(kk) qw]=KLDiv(h2,h0);
         [qq KLEntropy(kk)]=KLDiv(h0,h2);
         [KLDivergenceToData(kk) KLEntropyToData(kk)]=KLDiv(h0,h1);
     end
      %assignin('base', 'KLDivergence', KLDivergence);
      figure ('position', [00, 10, 650, 500])
       plot([1:OtherPara.TrajLength(Id_index)],KLDivergence,'.','markersize',20,'color','b'); hold on;
         plot([1:OtherPara.TrajLength(Id_index)],KLEntropy,'.','markersize',20,'color','r');
       ylabel('KL Divergence and entropy');
       xlabel('Time (h)');%ylabel('Single cells');
    %xticks([1:24:150]);xticklabels({'0','2','4','6','8','10','12'});
       xlabel('Time (h)');
       if OtherPara.Dataset==2 % 2 is p53 after Mdmx;
         xticks([1:20:122]);xticklabels({'0','10','20','30','40','50','60'});
         %xticks([1:24:130]);xticklabels({'0','2','4','6','8','10','12'});
    elseif OtherPara.Dataset==3 % 3 is Erk
         xticks([1:40:130]);xticklabels({'0','5','10','15'});
         %xticks([1:24:130]);xticklabels({'0','2','4','6','8','10','12'});
       else%  elseif OtherPara.Dataset==1
        xticks([1:24:130]);xticklabels({'0','2','4','6','8','10','12'});
    
        xlim([1 max(OtherPara.TimePointsToUse)])
        TimeUnit=round(OtherPara.TotalTimeLength/6);
        xTime=0:TimeUnit:5*TimeUnit;
        xTimePoints=round(xTime*OtherPara.TimePointsToUse/OtherPara.TotalTimeLength);    
        xTimePoints(1)=1;
        xticks(xTimePoints);xticklabels(cellstr(string(xTime)));%xticklabels({'0','2','4','6','8','10','12'});
    end
       title('Distance between sample and data');
       h=legend('KL Divergence','Entropy of distribution','Location','best');
       set(h,'FontSize',16);
        set(gca,'FontSize',20); 
        figurenamehmm=[foldername,'\TimeDepenKLDivergence_ToData_X',num2str(OtherPara.EnhancedSample),'_',num2str(ID(Id_index)),'States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize),'.jpg'];   
         %%saveas(gcf,figurenamehmm);
        
         HMM.KLDivergence{Id_index}=KLDivergence;HMM.KLEntropy{Id_index}=KLEntropy;
         
         
         figure ('position', [00, 10, 650, 500])
       plot([1:OtherPara.TrajLength(Id_index)],KLDivergence./KLEntropy,'.','markersize',20); hold on;
       ylabel('KL Divergence / entropy');xlabel('Time (h)');xlabel('Time (h)');
       if OtherPara.Dataset==2 % 2 is p53 after Mdmx;
         xticks([1:20:122]);xticklabels({'0','10','20','30','40','50','60'});
         %xticks([1:24:130]);xticklabels({'0','2','4','6','8','10','12'});
    elseif OtherPara.Dataset==3 % 3 is Erk
         xticks([1:40:130]);xticklabels({'0','5','10','15'});
         %xticks([1:24:130]);xticklabels({'0','2','4','6','8','10','12'});
       else%  elseif OtherPara.Dataset==1
        xticks([1:24:130]);xticklabels({'0','2','4','6','8','10','12'});
   
        xlim([1 max(OtherPara.TimePointsToUse)])
        TimeUnit=round(OtherPara.TotalTimeLength/6);
        xTime=0:TimeUnit:5*TimeUnit;
        xTimePoints=round(xTime*OtherPara.TimePointsToUse/OtherPara.TotalTimeLength);    
        xTimePoints(1)=1;
        xticks(xTimePoints);xticklabels(cellstr(string(xTime)));%xticklabels({'0','2','4','6','8','10','12'});
       end
    title('Distance between sample and data');
       ylim([0 1]);
        set(gca,'FontSize',20); 
        figurenamehmm=[foldername,'\TimeDepenKLRatio_ToData_X',num2str(OtherPara.EnhancedSample),'_',num2str(ID(Id_index)),'States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize),'.jpg'];   
         %%saveas(gcf,figurenamehmm);
         
         
         
         
%          figure ('position', [00, 10, 650, 500])
%        plot([1:OtherPara.TrajLength(Id_index)]*5/60,KLDivergenceToData,'.','markersize',20,'color','b'); hold on;
%          plot([1:OtherPara.TrajLength(Id_index)]*5/60,KLEntropyToData,'.','markersize',20,'color','r');
%        ylabel('KL Divergence and entropy');xlabel('Time (h)');xlabel('Time (h)');xticks([0:2:12]);% xticks([1:24:150]);xticklabels({'0','2','4','6','8','10','12'});%xticks([0:2:12]);xlim([0,130/12]);%to be consistent with heatmap%xlabel('Time points');
%        title('Distance between binned and data');
%        h=legend('KL Divergence','Entropy of distribution','Location','best');
%        set(h,'FontSize',16);
%         set(gca,'FontSize',20); 
%         figurenamehmm=[foldername,'\TimeDepenKLDivergence_BinTodata_',num2str(ID(Id_index)),'States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize),'.jpg'];   
%          %%saveas(gcf,figurenamehmm);
%         
%          
%          figure ('position', [00, 10, 650, 500])
%        plot([1:OtherPara.TrajLength(Id_index)]*5/60,KLDivergenceToData./KLEntropyToData,'.','markersize',20); hold on;
%        ylabel('KL Divergence / entropy');xlabel('Time (h)');xlabel('Time (h)');xticks([0:2:12]);% xticks([1:24:150]);xticklabels({'0','2','4','6','8','10','12'});%xticks([0:2:12]);xlim([0,130/12]);%to be consistent with heatmap%xlabel('Time points');
%         title('Distance between binned and data');
%        ylim([0 1]);
%         set(gca,'FontSize',20); 
%         figurenamehmm=[foldername,'\TimeDepenKLRatio_BinToData_',num2str(ID(Id_index)),'States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize),'.jpg'];   
%          %saveas(gcf,figurenamehmm);
         
       
     close all;
     continue;
 end
  if OtherPara.CalculateDistance==1 && OtherPara.TrajBasedStatistic==1
          
      
            
         [metrics,aux] = nfkbmetrics(Y{Id_index}');
         [aaa,indexData]=sortrows(metrics.peakfreq);
         indexData=flipud(indexData);
        
         
         if OtherPara.HMMFitToAll==1
             figurenamehmm=[foldername,'\Dist_X',num2str(OtherPara.EnhancedSample),'_',num2str(ID(Id_index)),'States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize),'.jpg'];
         else
             TrajSample_Smooth{1,Id_index}=SampledTraj{1,Id_index};
             OtherPara.TrajLength(Id_index)=size(SampledTraj{1,Id_index},2);
             figurenamehmm=[foldername,'\Dist_',num2str(ID(Id_index)),'States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize),'.jpg'];
         end
         
          [metrics2,aux] = nfkbmetrics(seqs{Id_index});
         [aaa,indexData2]=sortrows(metrics2.peakfreq);
         indexData2=flipud(indexData2);
         [metrics3,aux] = nfkbmetrics(TrajSample_Smooth{Id_index});
         [aaa,indexData3]=sortrows(metrics3.peakfreq);
         indexData3=flipud(indexData3);
       
      DistName={'euclidean','correlation','spearman'};
      
      % disp(TrajSample_Smooth);% disp(Id_index);disp(Y);
        IndexKSummary=[];
        for ll=1:size(TrajSample_Smooth{1,Id_index},1)
            if OtherPara.HMMFitToAll~=2
                [IndexK,DistK]=annquery(Y{Id_index},TrajSample_Smooth{1,Id_index}(ll,:)',1); %the first is data as reference
                IndexK2=ll;
            else
                IndexK=indexData(ll);
                IndexK2=indexData3(ll);
            end
      IndexKSummary(ll)=IndexK;
        %disp(IndexK);
         DistanceTrajs{Id_index}(1,ll)=pdist2(Y{Id_index}(:,IndexK)',TrajSample_Smooth{1,Id_index}(IndexK2,1:OtherPara.TrajLength(Id_index)),'euclidean');%It is already mean square
         DistanceTrajs{Id_index}(2,ll)=pdist2(Y{Id_index}(:,IndexK)',TrajSample_Smooth{1,Id_index}(IndexK2,1:OtherPara.TrajLength(Id_index)),'correlation');
         DistanceTrajs{Id_index}(3,ll)=pdist2(Y{Id_index}(:,IndexK)',TrajSample_Smooth{1,Id_index}(IndexK2,1:OtherPara.TrajLength(Id_index)),'spearman');
        end
        DistanceTrajs{Id_index}(1,:)=sqrt(DistanceTrajs{Id_index}(1,:).^2/size(Y{Id_index},1))/mean(mean(Y{Id_index})); % Normalize Euclidean to be CV-like
        figure ('position', [00, 10, 650, 500])
        for xx=1:1%3
        %subplot(1,3,xx)
         histogram(DistanceTrajs{Id_index}(xx,:),15,'FaceColor','b','EdgeColor','b');
         DistanceConditions(Id_index,xx)=mean(DistanceTrajs{Id_index}(xx,:),2);
        h=legend(num2str(round(DistanceConditions(Id_index,xx),2)),'Location','northeast');
        set(h,'FontSize',12);
        xlim([0 1.5])
        title(DistName(xx));
        set(gca,'FontSize',20);
        end
        figurenamehmm=[foldername,'\Dist_',num2str(OtherPara.EnhancedSample),'_',num2str(ID(Id_index)),'States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize),'.jpg'];        
         %saveas(gcf,figurenamehmm);
        
      
         
         
         %%
          DistName={'euclidean','correlation','spearman'};
        
        for ll=1:size(TrajSample_Smooth{1,Id_index},1)
          if OtherPara.HMMFitToAll~=2
              [IndexK,DistK]=annquery(seqs{Id_index}',TrajSample_Smooth{1,Id_index}(ll,:)',1); %the first is data as reference
              IndexK2=ll;
          else
              IndexK=indexData2(ll);
              IndexK2=indexData3(ll);
          end
         
         DistanceTrajs{Id_index}(1,ll)=pdist2(seqs{Id_index}(IndexK,:),TrajSample_Smooth{1,Id_index}(IndexK2,1:OtherPara.TrajLength(Id_index)),'euclidean');%It is already mean square
         DistanceTrajs{Id_index}(2,ll)=pdist2(seqs{Id_index}(IndexK,:),TrajSample_Smooth{1,Id_index}(IndexK2,1:OtherPara.TrajLength(Id_index)),'correlation');
         DistanceTrajs{Id_index}(3,ll)=pdist2(seqs{Id_index}(IndexK,:),TrajSample_Smooth{1,Id_index}(IndexK2,1:OtherPara.TrajLength(Id_index)),'spearman');
          
        end
        HMM.BiasRatio(Id_index)=length(unique(IndexKSummary))/length(IndexKSummary);
       
%         %assignin('base', 'test0', HMM.BiasRatio);
% %         %assignin('base', 'test0', sort(IndexKSummary));     
% %          %assignin('base', 'test1', seqs{Id_index}');
% %          %assignin('base', 'test2',TrajSample_Smooth{1,Id_index}');
% %          ddd

        DistanceTrajs{Id_index}(1,:)=sqrt(DistanceTrajs{Id_index}(1,:).^2/size(seqs{Id_index},2))/mean(mean(seqs{Id_index})); % Normalize Euclidean to be CV-like
        figure ('position', [00, 10, 650, 500])
        for xx=1:1%3
        %subplot(1,3,xx)
         histogram(DistanceTrajs{Id_index}(xx,:),15,'FaceColor','b','EdgeColor','b');
         DistanceConditions(Id_index,xx)=mean(DistanceTrajs{Id_index}(xx,:),2);
        h=legend(num2str(round(DistanceConditions(Id_index,xx),2)),'Location','northeast');
        set(h,'FontSize',12);
        xlim([0 1.5])
        title(DistName(xx));
        set(gca,'FontSize',20);
        end
        figurenamehmm=[foldername,'\Dist_ToBin_',num2str(ID(Id_index)),'States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize),'.jpg'];   
         %saveas(gcf,figurenamehmm);
         
         
         
        %% 
         %relative error
        
         for ll=1:size(TrajSample_Smooth{1,Id_index},1)
           %[IndexK,DistK]=annquery(seqs{Id_index}',TrajSample_Smooth{1,Id_index}(ll,:)',1); %the first is data as reference       
           if OtherPara.HMMFitToAll~=2
               [IndexK,DistK]=annquery(seqs{Id_index}',TrajSample_Smooth{1,Id_index}(ll,:)',1); %the first is data as reference
               IndexK2=ll;
           else
               IndexK=indexData2(ll);
               IndexK2=indexData3(ll);
           end
            %Distance2_Trajs{Id_index}(1,ll)=mean(sqrt((seqs{Id_index}(IndexK,:)-TrajSample_Smooth{1,Id_index}(IndexK2,1:OtherPara.TrajLength(Id_index))).^2)./seqs{Id_index}(IndexK,:));%relative error
            
            Distance2_Trajs{Id_index}(ll,:)=sqrt((seqs{Id_index}(IndexK,:)-TrajSample_Smooth{1,Id_index}(IndexK2,1:OtherPara.TrajLength(Id_index))).^2)./seqs{Id_index}(IndexK,:);%relative error
        end
        figure ('position', [00, 10, 650, 500])
        for xx=1:1%3
        %subplot(1,3,xx)
         histogram(mean(Distance2_Trajs{Id_index},2),15,'FaceColor','b','EdgeColor','b');
         DistanceConditions(Id_index,xx)=mean(mean(Distance2_Trajs{Id_index},2));
        h=legend(num2str(round(DistanceConditions(Id_index,xx),2)),'Location','northeast');
        set(h,'FontSize',12);
        xlim([0 1.5])
        title(DistName(xx));
        set(gca,'FontSize',20);
        end
        figurenamehmm=[foldername,'\RelError_ToBin_',num2str(ID(Id_index)),'States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize),'.jpg'];   
         %saveas(gcf,figurenamehmm);
         
         figure ('position', [00, 10, 650, 500])
        for xx=1:1%3
        %subplot(1,3,xx)
         %errorbar(1:size(Distance2_Trajs{Id_index},2),mean(Distance2_Trajs{Id_index},1),std(Distance2_Trajs{Id_index},1),'.','markersize',10);
         plot(1:size(Distance2_Trajs{Id_index},2),mean(Distance2_Trajs{Id_index},1),'.','markersize',10);
         %DistanceConditions(Id_index,xx)=mean(mean(Distance2_Trajs{Id_index},2));
        %h=legend(num2str(round(DistanceConditions(Id_index,xx),2)),'Location','northeast');
        %set(h,'FontSize',12);
        ylim([0 1]);ylabel('Relative error');
        title(DistName(xx));
        set(gca,'FontSize',20);
        end
        figurenamehmm=[foldername,'\RelErrorTimeDepen_ToBin_',num2str(ID(Id_index)),'States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize),'.jpg'];   
         %saveas(gcf,figurenamehmm);
       
  end
    close all;

end
disp('Calculate distance done!');


 
 close all;