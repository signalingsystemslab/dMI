function MostProbable(estTR,estE,OtherPara,seqs,Y,Y2,Nq,foldername,ID,HMM,SampledTraj)
Y=Y2;
Nq=length(OtherPara.binsize);
disp(ID);
logpseq=cell(1,Nq);MostProb=cell(1,Nq);

for Id_index=1:Nq%+1
     %seqs{Id_index}=seqs{Id_index}*OtherPara.Conversion(Id_index);
     for j=1:size(seqs{Id_index},1)
      likelystates = hmmviterbi(seqs{Id_index}(j,:),estTR{Id_index},estE{Id_index});
      MostProb{Id_index}(j)=sum(seqs{Id_index}(j,:)==likelystates)/size(seqs{Id_index},2);
    %PosteriorProb{Id_index}(j,:)=sum(PSTATES,1);
     end
     
     assignin('base', 'likelystates', likelystates);
     assignin('base', 'MostProb', MostProb);
%   figurenamehmm=[foldername,'\Data_FitToAll_',num2str(ID(1)),'States_',num2str(OtherPara.state(Id_index)),'bin_',num2str(OtherPara.binsize(Id_index)),'.jpg'];
%  saveas(gcf,figurenamehmm); 
 
dddd
end
disp('Calculate logP done!');


 
 close all;