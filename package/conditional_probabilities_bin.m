function [TMarix]= conditional_probabilities_bin(XX,OtherPara)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% GETINFO calculates mutual information for:
% 
% q         a set of input probabilities (sum(q)=1)
% Fcond     an [nxn] conditional probability matrix
% hRS       a vector of conditional entropies
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
% Get rows of conditional probability matrix - one "query", all "references"
% Sum across this row, individually weighting each point by appropriate
% weight.
%d=size(Fcond1_All,2)/size(Fcond1_All,1)+1;
Nq=length(XX);%number of conditions
TMarix=cell(Nq,1);
filefolder=pwd;
%edges=linspace(OtherPara.MinValue,OtherPara.MaxValue,OtherPara.binsize);
 edges=linspace(OtherPara.MinValue,OtherPara.MaxValue,OtherPara.binsize); % should not +1 for hist3
global ID

 %figure

for i=1:Nq
    
   
    [N,c]=hist3(XX{i}','Edges',{edges edges},'CDataMode','auto','FaceColor','interp'); %XX{i} is 2*traj: bivariate histogram plot of X(:,1) and X(:,2)
    N=N/sum(sum(N));% normalize the joint distribution

    %     colorbar
%     view(2)
%     xlabel('NFkB pre-value')
%     title(num2str(OtherPara.k))
%     ylabel('NFkB post-value')
%     set(gca,'FontSize',20);
    
%  figurename2=[filefolder,'/TransitionMatrix/',num2str(ID(i)),'_',num2str(OtherPara.k),'.jpg'];
%  saveas(gcf,figurename2); 

 TMarix{i}=N;

end

%dd
 
return;

% for kk=1:d-1
% %display(kk);
%     
% 
% % Fcond1=Fcond1_All(:,(kk-1)*Nq+1:kk*Nq);
% % Fcond2=Fcond2_All(:,(kk-1)*Nq+1:kk*Nq);
% 
%  disp(Fcond2{1,1});dd
% % disp(Fcond2);disp(size(Fcond2));
% 
% for s1=1:Nq
% %hRS(s1)=-sum(log2(Fcond2{s1,s1}(Fcond1{s1,s1}>eps)./Fcond1{s1,s1}(Fcond1{s1,s1}>eps)))/nnz((Fcond1{s1,s1}>eps));
% if kk==1
%     hRSTemp{s1}=-log2(Fcond2{s1,s1}(Fcond1{s1,s1}>eps)./Fcond1{s1,s1}(Fcond1{s1,s1}>eps))/nnz((Fcond1{s1,s1}>eps));
% else
%     hRSTemp{s1}=hRSTemp{s1}-log2(Fcond2{s1,s1}(Fcond1{s1,s1}>eps)./Fcond1{s1,s1}(Fcond1{s1,s1}>eps))/nnz((Fcond1{s1,s1}>eps));
% end
% if kk==d-1
%     hRS(s1)=sum(hRSTemp{s1});
% end
% end
% 
% end

