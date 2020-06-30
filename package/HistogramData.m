% Plot histogram for all the measured values, for selecting a uniform threshold

function HistogramData(Avector,condition,foldername)


h2=figure ('position', [00, 10, 1700, 1000]);
cc=linspecer(length(Avector));
for i=1:length(Avector)
    %subplot(4,5,i)
    hhh=histogram( Avector{i},20,'FaceAlpha',0.6,'FaceColor',cc(i,:)); hold on;
    count1(i)=length(Avector{i});
    count2(i)=length(find(Avector{i}>10));
end
    yLimits = get(gca,'YLim');
    text(4,yLimits(2)*0.9,[num2str((1-round(sum(count2)/sum(count1),4))*100),'%'],'FontSize',30);
    text(14,yLimits(2)*0.9,[num2str(round(sum(count2)/sum(count1),4)*100),'%'],'FontSize',30);
    
    ylim([0 yLimits(2)*1.1]);
    plot([10 10],[0, 1e10],'linewidth',2,'color','b');
    ax = gca;
    ax.YRuler.Exponent = 0;
h=legend([string(condition),'Cutoff'],'Location','bestoutside');
%title(h,'States number');
set(h,'FontSize',26,'fontWeight','normal');
    xlabel('Values of data');
  ylabel('Counts of data');
 xlim([0 20])
 title('Histogram of data points','FontWeight','Normal');
set(gca,'FontSize',30,'linewidth',2);

figurenamehmm=[foldername,'\HistoOfData.jpg'];
saveas(gcf,figurenamehmm); 


% h2=figure ('position', [-1800, 10, 2000, 1500]);
% for i=1:length(Avector)
%     subplot(4,5,i)
%     hhh=histogram( Avector{i},20); hold on;
%     ylim([0 max(hhh.Values)*1.1]);
%     plot([10 10],[0, 1e10]);
%     ax = gca;
%     ax.YRuler.Exponent = 0;
% 
%     xlabel('Values of data');
%   ylabel('Counts of data');
%  xlim([0 20])
%  title(string(condition{i}),'FontWeight','Normal');
%  set(gca,'FontSize',16);
% end
% figurenamehmm=[foldername,'\HistoOfData.jpg'];
% saveas(gcf,figurenamehmm); 
end