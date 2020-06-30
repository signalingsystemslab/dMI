%Function to calculate channel-capacity

function [MI_Full1 aSummaryTemp]= getInfoHMM_New(q,HMM,Nq1,OtherPara)
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
q=real(q);

if OtherPara.Jackknife==1
    TrajProbToUse=HMM.tempTrajProb;
else
    TrajProbToUse=HMM.TrajProb;
end

% if Nq1==6 %plot 6 conditions for IkB-mutant comparison, to make the first 6 conditions used here more informative
%    % disp('delete conditions......')
%     DeleteNum=1:2:20;
%     TrajProbToUse(DeleteNum,:)=[];
%     TrajProbToUse(:,DeleteNum)=[];
% end

%assignin('base', 'aa4', TrajProbToUse);
if OtherPara.CCSelect==1 %calculate CC for selected conditions, such as pairwise or 4 conditions etc
   % disp('delete conditions......')
    DeleteNum=1:min(size(HMM.TrajProb,1),30);
    %disp(size(HMM.TrajProb,1));dddd%1:26;
    assignin('base', 'aa2', DeleteNum);
    DeleteNum(OtherPara.CCSelectCondi)=[];
    TrajProbToUse(DeleteNum,:)=[];
    TrajProbToUse(:,DeleteNum)=[];
    Nq1=length(OtherPara.CCSelectCondi);
    %assignin('base', 'aa3', TrajProbToUse);
    %assignin('base', 'aa3', Nq1);
    
end

%assignin('base', 'aa', TrajProbToUse);




NumTime=size(TrajProbToUse{1,1},2);
% 

Nq=size(TrajProbToUse,2);%number of conditions
MI_Full=zeros(1,NumTime);

ConditionsAverageTrajProb=cell(1,Nq1);
CompensationFactor=0;
for ii=1:Nq1%Nq
    %Prefactor=[];
    %CompensationFactor=CompensationFactor+size(TrajProbToUse{1,ii},1);
   %sum(log(HMM.TrajEigen1{Id_index}(jjj,1:kk)+realmin))/kk;
    % ConditionsAverageTrajProb(ii,:)=mean(HMM.TrajEigen1{ii,ii},1);
    %TrajNum(ii)=size(TrajProbToUse{1,ii},1);
    for jj=1:Nq1
        if jj==1
            ConditionsAverageTrajProb{1,ii}(:,OtherPara.Point)=real(TrajProbToUse{jj,ii}(:,OtherPara.Point)*q(jj)); 
            %ConditionsAverageTrajEntropyFiltered(ii,:)=real(mean(HMM.TrajEntropyFiltered{jj,ii},1)*q(jj));     
        else
           % ConditionsAverageTrajProb{1,ii}(:,OtherPara.Point)=ConditionsAverageTrajProb{1,ii}(:,OtherPara.Point)+real(TrajProbToUse{jj,ii}(:,OtherPara.Point)*q(jj)); 
            %It seems that occasionally the above property breaks down for
            %some trajectory, thus remove those by making it to the
            %unconditional probability:
            Penalty=10;
            ConditionsAverageTrajProb{1,ii}(:,OtherPara.Point)=ConditionsAverageTrajProb{1,ii}(:,OtherPara.Point)+real(min(TrajProbToUse{jj,ii}(:,OtherPara.Point),TrajProbToUse{ii,ii}(:,OtherPara.Point)*Penalty)*q(jj)); 
            
            %ConditionsAverageTrajEntropyFiltered(ii,:)=ConditionsAverageTrajEntropyFiltered(ii,:)+real(mean(HMM.TrajEntropyFiltered{jj,ii},1)*q(jj)); 
    
        end
    end
    %HrTemp(ii,:)=real(mean(-(ConditionsAverageTrajProb{1,ii}.*log2(ConditionsAverageTrajProb{1,ii}+realmin))/size(ConditionsAverageTrajProb{1,ii},1)));
%    for kk=1:size(ConditionsAverageTrajProb{1,ii},1)
%     Prefactor(kk,:)=1./[1:1:size(ConditionsAverageTrajProb{1,ii},2)];
%    end
   %HrTemp(ii,:)=real(mean(-((Prefactor).*log2(ConditionsAverageTrajProb{1,ii})),1));
   %HrTemp(ii,OtherPara.Point)=(1/OtherPara.Point)*real(mean(-(log2(ConditionsAverageTrajProb{1,ii}(:,OtherPara.Point))),1));
   HrTemp(ii,OtherPara.Point)=real(mean(-(log2(ConditionsAverageTrajProb{1,ii}(:,OtherPara.Point))),1));%-OtherPara.Point*log2(OtherPara.binsize);
   %HrTemp(ii,OtherPara.Point)=real(mean(-(ConditionsAverageTrajProb{1,ii}(:,OtherPara.Point).*log2(ConditionsAverageTrajProb{1,ii}(:,OtherPara.Point))),1));%-OtherPara.Point*log2(OtherPara.binsize);
   
   %HrTemp(ii,:)=real(mean(-(log2(ConditionsAverageTrajProb{1,ii}+realmin))/size(ConditionsAverageTrajProb{1,ii},1)));
   % disp(indexTemp);
    %hRS(ii,:)=mean(-HMM.TrajLyapunov1{ii,ii},1);
    %hRS(ii,:)=real(mean(HMM.TrajEntropy{ii,ii},1));
    %Prefactor2(ii,:)=1:1:length(mean(HMM.TrajEntropy{ii,ii},1));
    %hRSTemp(ii,OtherPara.Point)=real(OtherPara.Point.*mean(HMM.TrajEntropy{ii,ii}(:,OtherPara.Point),1));
    %hRSTemp(ii,OtherPara.Point)=(1/OtherPara.Point)*real(mean(-log2(TrajProbToUse{ii,ii}(:,OtherPara.Point)),1));
    hRSTemp(ii,OtherPara.Point)=real(mean(-log2(TrajProbToUse{ii,ii}(:,OtherPara.Point)),1));%-OtherPara.Point*log2(OtherPara.binsize);
    %hRSTemp(ii,OtherPara.Point)=real(mean(-TrajProbToUse{ii,ii}(:,OtherPara.Point).*log2(TrajProbToUse{ii,ii}(:,OtherPara.Point)),1));%-OtherPara.Point*log2(OtherPara.binsize);
    
    %hRS_2(ii,:)=real(mean(HMM.TrajEntropyFiltered{ii,ii},1));
   % hRS(ii,:)=mean(-(TrajProbToUse{ii,ii}).*log2(TrajProbToUse{ii,ii}+realmin),1);
    %disp(mean(hRS(ii,:),2));
     
end

% %assignin('base', 'HrTemp', HrTemp); 
% %assignin('base', 'Prefactor', Prefactor); 
% %assignin('base', 'HMM1', TrajProbToUse); 

% if isfield(OtherPara,'CompensationFactor')
%     CompensationFactor=OtherPara.CompensationFactor;
% else
%     CompensationFactor=CompensationFactor/Nq1;
% end
% CompensationFactor=642.5000;
%disp('For our data, manually fix CompensationFactor in getInfoHMM');
%disp(CompensationFactor);
%HrTemp=HrTemp*CompensationFactor;
%dd
% %assignin('base', 'ConditionsAverageTrajProb', ConditionsAverageTrajProb);    
% dd
%disp(sum(ConditionsAverageTrajProb,2));
%Hr=0;%hRSTemp=cell(1,Nq);
for kk=OtherPara.Point%OtherPara.Point%1:NumTime %Do path average step by step...
%display(kk);



%Temp(:,kk)=-log2(ConditionsAverageTrajProb(:,kk)+realmin)/kk;

%Temp2=ConditionsAverageTrajProb(:,kk)'*q;
%disp(size(ConditionsAverageTrajProb)); disp(size(q));
% Temp2=ConditionsAverageTrajProb(:,kk)'*q;
% 
% Temp3=cell(1,Nq1);
% for ii=1:Nq1
% Temp3{ii}=TrajProbToUse{ii,ii}(:,kk)*q(ii);
% HrTemp(ii)=-sum(log2(Temp3{ii}+realmin))/kk;
% end

%Hr(kk)=-log2(ConditionsAverageTrajProb(:,kk)+realmin)'*q;
%Hr(kk)=real(-((HrTemp(:,kk)).*log2(HrTemp(:,kk)+realmin))'*q);
%Hr(kk)=real(-(1/CompensationFactor.*log2(ConditionsAverageTrajProb(:,kk)+realmin))'*q);
%Hr(kk)=real(-(1./TrajNum'.*log2(ConditionsAverageTrajProb(:,kk)+realmin))'*q);
%Hr(kk)=real(-sum((ConditionsAverageTrajProb(:,kk)).*log2(ConditionsAverageTrajProb(:,kk)+realmin),1));
Hr(kk)=real(HrTemp(:,kk)'*q);

%Hr_2(kk)=real(ConditionsAverageTrajEntropyFiltered(:,kk)'*q);

%disp(size(Temp2));disp(size(Temp2(Temp2>realmin)*log2(Temp2(Temp2>realmin))));
% if size(Temp2(Temp2>realmin),1)>0
%     Hr(kk)=-Temp2(Temp2>realmin)*log2(Temp2(Temp2>realmin));  
% else
%     Hr(kk)=0;  
% end




%Hr(kk)=-sum(log2(Temp2+realmin))/kk;
%Hr(kk)=sum(HrTemp);
%Hr(kk)=-log2(Temp4+realmin)/kk;


hRS(kk)=real(hRSTemp(:,kk)'*q);

MI_Full(kk)=Hr(kk)-hRS(kk);

%MI_Full_2(kk)=Hr_2(kk)-hRS_2(:,kk)'*q;

% disp(f_temp);dd
end
%disp(MI_Full);

% if min(MI_Full)<-2
%     dd
% end


 
 
%%
% Entropy = difference of H(R) and H(R|S), summed over the probabilities of each response. 
%f=real(sum(MI_Full));
f=real(MI_Full(OtherPara.Point));
MI_Full1=real(MI_Full(OtherPara.Point));
%MI_Full1=abs(real(MI_Full(OtherPara.Point)));
%MI_Full1=real(sum(MI_Full));

% f=abs(real(MI_Full(OtherPara.Point)));


% disp(q); disp(HrTemp);
aSummaryTemp.q=q(:,end)';
aSummaryTemp.a0=HrTemp;
aSummaryTemp.a1=Hr;
aSummaryTemp.a2=hRSTemp;
aSummaryTemp.a3=hRS;
% %assignin('base', 'a0', HrTemp);
% %assignin('base', 'a1', Hr);
% %assignin('base', 'a2', hRSTemp);
% %assignin('base', 'a3', hRS);
% if OtherPara.Point==95 || isnan(q(1))
% a0=ConditionsAverageTrajProb;
% a1=HrTemp;%-((ConditionsAverageTrajProb(:,:)).*log2(ConditionsAverageTrajProb(:,:)+realmin));
% a2=hRS;
% a3=-sum((log2(HrTemp(:,:)+realmin)),2);
% a4=hRSTemp;
% 
% %assignin('base', 'a0', a0);
% %assignin('base', 'a1', a1);
% %assignin('base', 'a2', a2);
% %assignin('base', 'a1sum', sum(a1,2));
% %assignin('base', 'a2sum', sum(a2,2));
% %assignin('base', 'a3', a3);
% %assignin('base', 'a4', a4);
% %assignin('base', 'a5', q);      
%            % ddd
% end
        
 

%if max(isnan(q))>0
    
%     %assignin('base', 'a0', HrTemp);
% %assignin('base', 'a1', Hr);
% %assignin('base', 'a2', hRSTemp);
% %assignin('base', 'a3', hRS);
% %assignin('base', 'a4', q);
% %assignin('base', 'a5', ConditionsAverageTrajProb);
%dddd
%end


 %disp(sum(Hr)); disp(sum(hRS(:,:)'*q)); 
%dd
% f=sum(Hr)-hRS*q;
%f_total=sum(f*timeunit);
