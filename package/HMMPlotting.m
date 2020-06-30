%Plot the trajectory entropy for each condition and time-dependent channel
%capacity, etc. More fine-tuning figures are plotted by separate scripts as
%specified in README.txt.

function HMMPlotting(HMM,seqs,Nq,ID,OtherPara,foldername,SampledTraj)
if OtherPara.HMMFitToAll==1
if OtherPara.TrajBasedStatistic~=2

 figureSize1=round(sqrt(Nq));
 if figureSize1~=1
 figureSize2=figureSize1+1;
 Plotsize1=1200;
 Plotsize2=800;
 else
 figureSize2=figureSize1+1;
 Plotsize1=400;
 Plotsize2=300;
 end
 
     h1=figure ('position', [00, 10, Plotsize1, Plotsize2]);
for Id_index=1:Nq%+1
subplot(figureSize2,figureSize1,Id_index)
if size(HMM.TrajProb{Id_index,Id_index},1)>1
errorbar(1:size(HMM.TrajProb{Id_index,Id_index},2),mean(HMM.TrajProb{Id_index,Id_index},1),std(HMM.TrajProb{Id_index,Id_index},1),'.','markersize',10); hold on;
else
plot(1:size(HMM.TrajProb{Id_index,Id_index},2),mean(HMM.TrajProb{Id_index,Id_index},1),'.','markersize',10); hold on;
end%plot(1:size(seqs{Nq},2),mean(HMM.TrajProb{Id_index,Id_index},1),'.','markersize',10); hold on;
%plot(1:size(seqs{Nq},2),sum(HMM.TrajEntropy{Id_index,Id_index},1),'.','markersize',10); hold on;
if Id_index<=Nq
    title(num2str(ID(Id_index)));
end
ylim([0 1])
set(gca,'FontSize',16,'linewidth',2);
end

% if OtherPara.Dataset==2 % 2 is p53 after Mdmx;
%          xticks([1:20:122]);xticklabels({'0','10','20','30','40','50','60'});
%          %xticks([1:24:130]);xticklabels({'0','2','4','6','8','10','12'});
%     elseif OtherPara.Dataset==3 % 3 is Erk
%          xticks([1:40:130]);xticklabels({'0','5','10','15'});
%          %xticks([1:24:130]);xticklabels({'0','2','4','6','8','10','12'});
%       elseif OtherPara.Dataset==1
%         xticks([1:24:130]);xticklabels({'0','2','4','6','8','10','12'});
%     else
%         xlim([1 max(OtherPara.TimePointsToUse)])
%         TimeUnit=round(OtherPara.TotalTimeLength/6);
%         xTime=0:TimeUnit:5*TimeUnit;
%         xTimePoints=round(xTime*OtherPara.TimePointsToUse/OtherPara.TotalTimeLength);    
%         xTimePoints(1)=1;
%         xticks(xTimePoints);xticklabels(cellstr(string(xTime)));%xticklabels({'0','2','4','6','8','10','12'});
% end
%     xlabel('Time (h)');
figurenamehmm=[foldername,'\Traj_Prob_States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize),'.jpg'];
%saveas(gcf,figurenamehmm);
% 
% h1=figure ('position', [00, 10, Plotsize1, Plotsize2]);
% for Id_index=1:Nq%+1
% subplot(figureSize2,figureSize1,Id_index)
% %disp(size(seqs{Nq},2));
% %disp(size(HMM.TrajEntropy{Id_index,Id_index}));
% %errorbar(1:size(seqs{Nq},2),mean(HMM.TrajEntropy{Id_index,Id_index},1),std(HMM.TrajEntropy{Id_index,Id_index},1),'.','markersize',10); hold on;
% %plot(1:size(HMM.TrajEntropy{Id_index,Id_index},2),mean(HMM.TrajEntropy{Id_index,Id_index},1),'.-','markersize',10); hold on;
% if size(HMM.TrajProb{Id_index,Id_index},1)>1
%     errorbar(1:size(HMM.TrajEntropy{Id_index,Id_index},2),mean(HMM.TrajEntropy{Id_index,Id_index},1),std(HMM.TrajEntropy{Id_index,Id_index},1),'.','markersize',10); hold on;
% else
%     plot(1:size(HMM.TrajEntropy{Id_index,Id_index},2),mean(HMM.TrajEntropy{Id_index,Id_index},1),'.','markersize',10); hold on;
% 
% end%plot(1:size(seqs{Nq},2),sum(HMM.TrajEntropy{Id_index,Id_index},1),'.','markersize',10); hold on;
% if Id_index<=Nq
%     title(num2str(ID(Id_index)));
% end
% ylim([0 5])
% set(gca,'FontSize',16,'linewidth',2);
% end
% figurenamehmm=[foldername,'\Traj_EntropyRate_States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize),'.jpg'];
% %%saveas(gcf,figurenamehmm);
% 
% 
% %  h1=figure ('position', [00, 10, Plotsize1, Plotsize2]);
% % for Id_index=1:Nq%+1
% % subplot(figureSize,figureSize-1,Id_index)
% % errorbar(1:size(HMM.TrajProb_Point{Id_index,Id_index},2),mean(HMM.TrajProb_Point{Id_index,Id_index},1),std(HMM.TrajProb_Point{Id_index,Id_index},1),'.','markersize',10); hold on;
% % %plot(1:size(seqs{Nq},2),sum(HMM.TrajEntropy{Id_index,Id_index},1),'.','markersize',10); hold on;
% % if Id_index<=Nq
% %     title(num2str(ID(Id_index)));
% % end
% % %ylim([0 0.5])
% % set(gca,'FontSize',16,'linewidth',2);
% % end
% % figurenamehmm=[foldername,'\TrajPoint_Prob_States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize),'.jpg'];
% % %%saveas(gcf,figurenamehmm);
% % 
% %  h1=figure ('position', [00, 10, Plotsize1, Plotsize2]);
% % for Id_index=1:Nq%+1
% % subplot(figureSize,figureSize-1,Id_index)
% % errorbar(1:size(HMM.TrajEntropy_Point{Id_index,Id_index},2),mean(HMM.TrajEntropy_Point{Id_index,Id_index},1),std(HMM.TrajEntropy_Point{Id_index,Id_index},1),'.','markersize',10); hold on;
% % %plot(1:size(seqs{Nq},2),sum(HMM.TrajEntropy{Id_index,Id_index},1),'.','markersize',10); hold on;
% % if Id_index<=Nq
% %     title(num2str(ID(Id_index)));
% % end
% % %ylim([0 0.5])
% % set(gca,'FontSize',16,'linewidth',2);
% % end
% % figurenamehmm=[foldername,'\TrajPoint_EntropyRate_States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize),'.jpg'];
% % %%saveas(gcf,figurenamehmm);
% 
% 
% 
% 
%  h1=figure ('position', [00, 10, Plotsize1, Plotsize2]);
% for Id_index=1:Nq%+1
% subplot(figureSize2,figureSize1,Id_index)
% %errorbar(1:size(seqs{Nq},2),mean(cumsum(HMM.TrajEntropy{Id_index,Id_index},2),1),std(cumsum(HMM.TrajEntropy{Id_index,Id_index},2),1),'.','markersize',10); hold on;
% plot(1:size(HMM.TrajEntropy{Id_index,Id_index},2),mean(cumsum(HMM.TrajEntropy{Id_index,Id_index},2),1),'.','markersize',10); hold on;
% if Id_index<=Nq
%     title(num2str(ID(Id_index)));
% end
% %ylim([0 0.5])
% set(gca,'FontSize',16,'linewidth',2);
% end
% figurenamehmm=[foldername,'\Traj_EntropyAccu_States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize),'.jpg'];
% %saveas(gcf,figurenamehmm);
% 
% 
% 
% %  h1=figure ('position', [00, 10, Plotsize1, Plotsize2]);
% % for Id_index=1:Nq%+1
% % subplot(figureSize2,figureSize1,Id_index)
% % errorbar(1:size(HMM.TrajEntropyFiltered{Id_index,Id_index},2),mean(HMM.TrajEntropyFiltered{Id_index,Id_index},1),std(HMM.TrajEntropyFiltered{Id_index,Id_index},1),'.','markersize',10); hold on;
% % %plot(1:size(seqs{Nq},2),cumsum(sum(HMM.TrajEntropy{Id_index,Id_index},1)),'.','markersize',10); hold on;
% % if Id_index<=Nq
% %     title(num2str(ID(Id_index)));
% % end
% % %ylim([0 0.5])
% % set(gca,'FontSize',16,'linewidth',2);
% % end
% % figurenamehmm=[foldername,'\TrajHiddenStates_Entropy_States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize),'.jpg'];
% % %saveas(gcf,figurenamehmm);
% 
% %  h1=figure ('position', [00, 10, Plotsize1, Plotsize2]);
% % for Id_index=1:Nq%+1
% % subplot(figureSize2,figureSize1,Id_index)
% % %plot(1:size(seqs{Nq},2),mean(HMM.TrajEntropyFiltered{Id_index,Id_index},1),std(HMM.TrajEntropyFiltered{Id_index,Id_index},1),'.','markersize',10); hold on;
% % plot(1:size(HMM.TrajEntropyFiltered{Id_index,Id_index},2),mean(cumsum(HMM.TrajEntropyFiltered{Id_index,Id_index},2),1),'.','markersize',10); hold on;
% % if Id_index<=Nq
% %     title(num2str(ID(Id_index)));
% % end
% % %ylim([0 0.5])
% % set(gca,'FontSize',16,'linewidth',2);
% % end
% % figurenamehmm=[foldername,'\TrajHiddenStates_EntropyAccu_States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize),'.jpg'];
% % %saveas(gcf,figurenamehmm);
% 
% 
% 
% h2=figure ('position', [00, 10, Plotsize1, Plotsize2]);
% for Id_index=1:Nq%+1
% subplot(figureSize2,figureSize1,Id_index)
% plot(1:size(HMM.TrajMemory{Id_index},2),mean(HMM.TrajMemory{Id_index},1),'.','markersize',10); hold on;
% %errorbar(1:size(HMM.TrajMemory{Id_index},2),mean(HMM.TrajMemory{Id_index},1),std(HMM.TrajMemory{Id_index},1),'.','markersize',10); hold on;
% %errorbar(1:size(seqs{Nq},2),mean(cumsum(HMM.TrajMemory{Id_index},2),1),std(cumsum(HMM.TrajMemory{Id_index},2),1),'.','markersize',10); hold on;
% if Id_index<=Nq
%     title(num2str(ID(Id_index)));
% end
% ylim([0 4])
% set(gca,'FontSize',16,'linewidth',2);
% end
% figurenamehmm=[foldername,'\Matrices_Memory_States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize),'.jpg'];
% %figurenamehmm=[foldername,'\Matrices_MemoryAccumu_States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize),'.jpg'];
% %saveas(gcf,figurenamehmm);
% 
%         
%     
% h2=figure ('position', [00, 10, Plotsize1, Plotsize2]);
% for Id_index=1:Nq%+1
% subplot(figureSize2,figureSize1,Id_index)
% if size(HMM.TrajProb{Id_index,Id_index},1)>1
% errorbar(1:size(HMM.TrajLyapunov1{Id_index},2),mean(HMM.TrajLyapunov1{Id_index},1),std(HMM.TrajLyapunov1{Id_index},1),'.','markersize',10); hold on;
% else
% plot(1:size(HMM.TrajLyapunov1{Id_index},2),mean(HMM.TrajLyapunov1{Id_index},1),'.','markersize',10); hold on;
% end
% if Id_index<=Nq
%     title(num2str(ID(Id_index)));
% end
% ylim([0 5])
% set(gca,'FontSize',16,'linewidth',2);
% end
% figurenamehmm=[foldername,'\Matrices_EntropyRate_States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize),'.jpg'];
% %saveas(gcf,figurenamehmm);


%HMM.MI_Full=abs(HMM.MI_Full);

h2=figure ('position', [00, 10, 400, 300]);
plot(1:size(HMM.MI_Full,2),max(HMM.MI_Full,0),'.','markersize',10); hold on;
ylim([0 5])
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
xlabel('Time (h)');
ylabel('Maximum MI (bits)');
set(gca,'FontSize',16,'linewidth',2);

if OtherPara.CCSelect==1
kk=1;Nowexist=1;
while Nowexist~=0
figurenamehmm=[foldername,'\CCofCondition_',num2str(kk),'_',num2str(Nq),'States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize),'.jpg'];
if exist(figurenamehmm)
    kk=kk+1;
else
    Nowexist=0;
end
end

else
figurenamehmm=[foldername,'\Channel_Condi_',num2str(Nq),'States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize),'.jpg'];
end

%saveas(gcf,figurenamehmm);
% h2=figure ('position', [00, 10, 400, 300]);
% plot(1:size(HMM.MI_Full,2),cumsum(HMM.MI_Full),'.','markersize',10); hold on;
% %ylim([0 2]);
% xlabel('Time points');
% ylabel('Accumulated channel capacity (bits)');
% set(gca,'FontSize',16,'linewidth',2);
% 
% figurenamehmm=[foldername,'\Matrices_ChannelAccumu_Condi_',num2str(Nq),'States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize),'.jpg'];
% %saveas(gcf,figurenamehmm);

MI.X=1:size(seqs{Nq},2);
MI.Y=HMM.MI_Full;


else
    MI=0;



end

return;
%% %model based static statistics
 if OtherPara.ModelBasedStatistic==1
     c = categorical(OtherPara.condition2(1:min(Nq,length(OtherPara.condition2))));
     

%      
% figure ('position', [00, 10, 1000, 400])
% for k=1:2
% subplot(1,2,k)
% c = categorical(OtherPara.condition2);
% if k==1
% barh(c,HMM.Memory(1:end-1))
% elseif k==2
% barh(c,HMM.Memory_2(1:end-1))
% end
% title('Memory decay rate');
% %ylim([0 0.5])
% set(gca,'FontSize',16,'linewidth',2);
% end
% figurenamehmm=[foldername,'\Memory_CondBased_States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize),'.jpg'];
% %saveas(gcf,figurenamehmm);
        
% figure ('position', [00, 10, 1000, 400])
% for k=1:2
% subplot(1,2,k)
% if k==1
% barh(c,HMM.Entropy(1:end))
% elseif k==2
% barh(c,HMM.Entropy_2(1:end))
% end
% title('Entropy');
% %ylim([0 0.5])
% set(gca,'FontSize',16,'linewidth',2);
% end
% figurenamehmm=[foldername,'\Entropy_CondBased_States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize),'.jpg'];
% %saveas(gcf,figurenamehmm);

figure ('position', [00, 10, 700, 800])
barh(c,HMM.Entropy(1:Nq))
title('Entropy');
%ylim([0 0.5])
set(gca,'FontSize',16,'linewidth',2);

figurenamehmm=[foldername,'\Entropy_CondBased_States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize),'.jpg'];
%saveas(gcf,figurenamehmm);


figure ('position', [00, 10, 700, 800])
barh(c,HMM.EntropyProduction(1:Nq))
title('Entropy Production');
%ylim([0 0.5])
set(gca,'FontSize',16,'linewidth',2);

figurenamehmm=[foldername,'\EntropyProd_CondBased_States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize),'.jpg'];
%saveas(gcf,figurenamehmm);


figure ('position', [00, 10, 700, 800])
TimeWindow=20;
for Id_index2=1:Nq
HMM.MemoryDecayCondition(Id_index2)=mean(mean(HMM.TrajMemory{Id_index2}(:,1:TimeWindow),1));
end
% %assignin('base', 'c1', c);
% %assignin('base', 'c2', HMM.MemoryDecayCondition);
barh(c,HMM.MemoryDecayCondition(1:Nq))
title('Average memory decay rate');
%ylim([0 0.5])
set(gca,'FontSize',16,'linewidth',2);

figurenamehmm=[foldername,'\MemoryDecay_Cond_States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize),'.jpg'];
%saveas(gcf,figurenamehmm);
% 
% bar.c=c;
% bar.EntropyProduction=HMM.EntropyProduction;
% bar.EntropyProduction=HMM.MemoryDecayCondition;
% 
%  else
%      bar=0;
 

     
     bar.c=c;
bar.EntropyProduction=HMM.EntropyProduction;
bar.MemoryDecayCondition=HMM.MemoryDecayCondition;
filename2=[foldername,'\Condi_',num2str(Nq),'States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize),'.mat'];
save(filename2,'MI','bar');
 end

%%
else
    %plot heat map
    TrajSample_Smooth=cell(1,Nq);
for Id_index=1:Nq

    if min(abs(ID(Id_index)-[548,546,778,780,720,610,779]))==0
        YUpLimit=5;
    elseif min(abs(ID(Id_index)-[650,566]))==0
        YUpLimit=4;
    elseif min(abs(ID(Id_index)-[756]))==0
        YUpLimit=7;
    else
        YUpLimit=OtherPara.MaxValue;
    end
if OtherPara.Dataset~=1 % 2 is p53 after Mdmx; 
    YUpLimit=OtherPara.MaxValue;% for p53 data
    end
opt= maketicks(1:size(SampledTraj{1,Id_index},2),[0, YUpLimit],0); opt.Name = 'ddd';
HM= figure ('position', [00, 10, 600, 850]);%('name','650');    

sampleSize=size(SampledTraj{1,Id_index},1);


for jjj=1:sampleSize
TrajSample_Smooth{1,Id_index}(jjj,:) = max(smooth(1:size(SampledTraj{1,Id_index},2),...
    SampledTraj{1,Id_index}(jjj,:),0.07,'loess'),0);
end
% %assignin('base', 'TrajSample_Smooth', TrajSample_Smooth); 
%[tblB,index]=sortrows(TrajSample_Smooth{1,Id_index},[1:1:120]);
%s=surf(TrajSample_Smooth{1,Id_index});
[metrics,aux] = nfkbmetrics(TrajSample_Smooth{1,Id_index});
[aaa,index]=sortrows(metrics.peakfreq);
tblB=TrajSample_Smooth{1,Id_index}(flipud(index),:);
    
h=colormapStack(tblB,[],opt, HM);
 hh = colorbar;  
    s.EdgeColor = 'none';view(2);
    xlabel('Time (h)');ylabel('Single cells');
   % zlim([0 4])
    if OtherPara.Dataset==2 % 2 is p53 after Mdmx;
         xticks([0:20:120]);xticklabels({'0','10','20','30','40','50','60'});
     else
    xticks([1:24:130]);xticklabels({'0','2','4','6','8','10','12'});
    end
    set(gca,'FontSize',40);
     figurenamehmm=[foldername,'\HMM_',num2str(ID(Id_index)),'States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize),'.jpg'];
 %saveas(gcf,figurenamehmm); 



 
 %DistanceTrajs=cell(1,Nq);
    
  

% if OtherPara.CalculateDistance==1
%     [metrics,aux] = nfkbmetrics(Y{Id_index}');
%          [aaa,indexData]=sortrows(metrics.peakfreq);
%          indexData=flipud(indexData);
%          
%       DistName={'euclidean','correlation','spearman'};
%         for ll=1:size(TrajSample_Smooth{1,Id_index},1)
%          [IndexK,DistK]=annquery(Y{Id_index},TrajSample_Smooth{1,Id_index}(ll,:)',1); %the first is data as reference
%          
%          DistK=indexData(ll);
%          DistanceTrajs{Id_index}(1,ll)=pdist2(Y{Id_index}(:,IndexK)',TrajSample_Smooth{1,Id_index}(ll,:),'euclidean');
%          DistanceTrajs{Id_index}(2,ll)=pdist2(Y{Id_index}(:,IndexK)',TrajSample_Smooth{1,Id_index}(ll,:),'correlation');
%          DistanceTrajs{Id_index}(3,ll)=pdist2(Y{Id_index}(:,IndexK)',TrajSample_Smooth{1,Id_index}(ll,:),'spearman');
%         end
%         DistanceTrajs{Id_index}(1,:)=DistanceTrajs{Id_index}(1,:)/mean(mean(Y{Id_index}))/sqrt(size(Y{Id_index},1)); % Normalize Euclidean to be CV-like
%         figure ('position', [00, 10, 900, 300])
%         for xx=1:3
%          subplot(1,3,xx)
%          histogram(DistanceTrajs{Id_index}(xx,:),15,'FaceColor','b','EdgeColor','b');
%          DistanceConditions(Id_index,xx)=mean(DistanceTrajs{Id_index}(xx,:),2);
%         h=legend(num2str(round(DistanceConditions(Id_index,xx),2)),'Location','northeast');
%         set(h,'FontSize',12);
%         title(DistName(xx));
%         set(gca,'FontSize',12);
%         end
%          figurenamehmm=[foldername,'\Dist_FitToAll_X',num2str(OtherPara.EnhancedSample),'_',num2str(ID(Id_index)),'States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize),'.jpg'];
%         %saveas(gcf,figurenamehmm);
% end
%     
end
%save(filename,'DistanceConditions','ID','DistanceTotal', 'SampledTraj','estTR','estE','estTR_sc','estE_sc','OtherPara');





end

