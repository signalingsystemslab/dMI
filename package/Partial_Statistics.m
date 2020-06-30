function [DistanceConditions,HMM,seqs]=Partial_Statistics(estTR,estE,OtherPara,Y,YTest,Nq,foldername,ID,HMM,SampledTraj)

Nq=length(OtherPara.binsize);
seqs=cell(1,Nq);seqs2=cell(1,Nq);
DistanceConditions=0;
    TrajSample=cell(1,Nq);
    TrajSample_Smooth=cell(1,Nq);
    TrajSample_HiddenStates=cell(1,Nq);
    estTRsort=cell(1,Nq);
    estEsort=cell(1,Nq);
    DistanceTrajs=cell(1,Nq);Distance2_Trajs=cell(1,Nq);
    
    
    if OtherPara.TrajBasedStatistic==1
    HMM.TrajProb=cell(Nq,Nq);HMM.TrajProb_Point=cell(Nq,Nq);HMM.TrajProb2=cell(Nq,Nq);
    HMM.TrajEntropy=cell(Nq,Nq);HMM.TrajEntropy_Point=cell(Nq,Nq);
    HMM.TrajMemory=cell(1,Nq);
    HMM.TrajEigen1=cell(1,Nq);HMM.TrajEigen2=cell(1,Nq);
    HMM.TrajLyapunov1=cell(Nq,Nq);HMM.TrajLyapunov2=cell(Nq,Nq);
    HMM.TrajEntropyFiltered=cell(Nq,Nq);
    end


IniDist=cell(Nq,1);
count1=0;count2=0;
for Id_index=1:Nq%+1
disp(Id_index);
        
seqSample_test=[];
 for jjj=1:OtherPara.PartialNum
       % disp(jjj);
    [seqSample2,states] = hmmgenerate(size(YTest{1},1),estTR{1,Id_index},estE{1,Id_index});
 
    seqSample_test(jjj,:)=seqSample2;%*OtherPara.Conversion(Id_index);
 end
 for type=1:1 % 1 is get Trajectory probability of test data, 2 is get Trajectory probability of sampled sets
    
 if type==1
    
    %Calculate model/condition based statistics
if OtherPara.HMMFitToAll==1 && OtherPara.CalculateStatistics==1 && OtherPara.TrajBasedStatistic==1
    %disp(Id_index);
    HMM.on=0;
    if OtherPara.ModelBasedStatistic==1
        [V,D] = eig(estTR{1,Id_index}');
        e= real(eig(estTR{1,Id_index}'));%-eye(size(estTR{1,Id_index},1)));
        [M,I]=min(abs(e-1));
        X=abs(real(V(:,I)));X=X/sum(X);
        %disp(X);
        if OtherPara.Dataset~=0
        X=ones(1,length(X));X=X/sum(X);% If using unifrom distribution
        end
        
        %X=linsolve(estTR{1,Id_index}'-eye(size(estTR{1,Id_index},1)),zeros(size(estTR{1,Id_index},1),1)); X=X/sum(X);disp(X);
       % X=rand(size(estTR{1,Id_index},1),1);X=X/sum(X);% Use random distribution 
       % X=ones(size(estTR{1,Id_index},1),1)/size(estTR{1,Id_index},1);% Use uniform distribution as steady state is all 0s.
        % Transition matrix only
        %MatrixTemp0=reshape(estTR{1,Id_index},[],1);
       % HMM.Entropy(Id_index)=sum(X'*estTR{1,Id_index}*logm(estTR{1,Id_index}));
        % Transition matrix .* E*E^tau
        MatrixTemp1=estTR{1,Id_index};
        %MatrixTemp2=estTR{1,Id_index}.*(estE{1,Id_index}*estE{1,Id_index}');
        %HMM.Entropy_2(Id_index)=sum(X'*MatrixTemp*logm(MatrixTemp));
        
%         HMM.Entropy(Id_index)=0;HMM.Entropy_2(Id_index)=0;
%         HMM.EntropyProduction(Id_index)=0;HMM.EntropyProduction_2(Id_index)=0;HMM.EntropyProduction00(Id_index)=0;
% %          assignin('base', 'MatrixTemp1', MatrixTemp1);     
% %          assignin('base', 'X', X);
% %dd
% %         EntropyProductionTest1=NaN(size(MatrixTemp1,1),size(MatrixTemp1,2));
% %         EntropyProductionTest2=NaN(size(MatrixTemp1,1),size(MatrixTemp1,2));
%          EntropyProductionTest=NaN(size(MatrixTemp1,1),size(MatrixTemp1,2));
%         for ii=1:size(MatrixTemp1,1)
%         for jj=1:size(MatrixTemp1,2)
%             if MatrixTemp1(ii,jj)>0
%                 HMM.Entropy(Id_index)=HMM.Entropy(Id_index)-X(ii)*MatrixTemp1(ii,jj)*log2(MatrixTemp1(ii,jj)+realmin);
%                 if MatrixTemp1(jj,ii)>0
%                     
% %                 HMM.EntropyProduction(Id_index)=HMM.EntropyProduction(Id_index)+X(ii)*MatrixTemp1(ii,jj)...
% %                     *log2(X(ii)*MatrixTemp1(ii,jj)/(X(jj)*MatrixTemp1(jj,ii)+realmin)+realmin);
%                 EntropyProductionTest(ii,jj)=MatrixTemp1(ii,jj)/MatrixTemp1(jj,ii);
%                 if X(ii)*MatrixTemp1(ii,jj)*log2(X(ii)*MatrixTemp1(ii,jj)/(X(jj)*MatrixTemp1(jj,ii)))<Inf
%                 HMM.EntropyProduction(Id_index)=HMM.EntropyProduction(Id_index)+X(ii)*MatrixTemp1(ii,jj)...
%                     *log2(X(ii)*MatrixTemp1(ii,jj)/(X(jj)*MatrixTemp1(jj,ii)));
%                 end
% %                 EntropyProductionTest1(ii,jj)=MatrixTemp1(ii,jj);
% %                 EntropyProductionTest2(ii,jj)=MatrixTemp1(jj,ii);
%                  
%                 
%                 
%                 
%                 HMM.EntropyProduction00(Id_index)=HMM.EntropyProduction00(Id_index)+X(ii)*MatrixTemp1(ii,jj)...
%                     *log2(MatrixTemp1(ii,jj)/(MatrixTemp1(jj,ii)+realmin)+realmin);
%                 end
%             end
% %             if MatrixTemp2(ii,jj)>0
% %                 HMM.Entropy_2(Id_index)=HMM.Entropy_2(Id_index)-X(ii)*MatrixTemp2(ii,jj)*log2(MatrixTemp2(ii,jj)+realmin);
% %                 if MatrixTemp2(jj,ii)>0
% %                 HMM.EntropyProduction_2(Id_index)=HMM.EntropyProduction_2(Id_index)+X(ii)*MatrixTemp2(ii,jj)...
% %                     *log2(X(ii)*MatrixTemp2(ii,jj)/(X(jj)*MatrixTemp2(jj,ii)+realmin)+realmin);
% %                 end
% %             end
%         end
%         end
%         assignin('base', 'EntropyProductionTest1',  X);
%         assignin('base', 'EntropyProductionTest2',  MatrixTemp1);
%          assignin('base', 'EntropyProductionTest3',  HMM.EntropyProduction);
%           assignin('base', 'EntropyProductionTest',  EntropyProductionTest);ddd
%         assignin('base', 'EntropyProductionTest1', EntropyProductionTest1);
%         assignin('base', 'EntropyProductionTest2', EntropyProductionTest2);
        
%         [Q,R] = qr(MatrixTemp1');e=diag(R);%e(e==0)=[]; 
%          e=sort(abs(e),'descend');
%          %e=sort(e,'descend');
%          %disp(e);
%        % e = sort( real(eig(MatrixTemp1')),'descend');
%         HMM.Eigen1(Id_index)=e(1); HMM.Eigen2(Id_index)=e(2);
%         HMM.Memory(Id_index)=abs(e(2)-e(1));
%         [Q,R] = qr(MatrixTemp2');e=diag(R);%e(e==0)=[]; 
%          e=sort(abs(e),'descend');
%          %  e=sort(e,'descend');
%         %e = sort( real(eig(MatrixTemp2')),'descend');
%         HMM.Eigen1_2(Id_index)=e(1); HMM.Eigen2_2(Id_index)=e(2);
%         HMM.Memory_2(Id_index)=abs(e(2)-e(1));      
    end
    
  %  continue;
    %%
     %Calculate traj based statistics
    if OtherPara.TrajBasedStatistic==1
        
         DeleteRowsNaN=[];DeleteRowsZero=[];
            
            
        if Id_index==Nq+1
         seqsUse=cat(1,seqs{1:Nq});
         YUse=cat(2,YTest{1:Nq});
        else
         
         YUse=YTest{1};
        end
        %assignin('base', 'YUse', YUse);
        for jj=1:size(YUse,2)
        [N,edgesNew,seqs{Id_index}(jj,:)] = histcounts(YUse(:,jj)',OtherPara.edges{Id_index});
         end

         seqsUse=seqs{Id_index};
         %assignin('base', 'seqsUse', seqsUse);
         
         [N,edgesNew,bin] = histcounts(YUse(1,:),OtherPara.edges{Id_index});
         
        %edges=linspace(OtherPara.MinValue,OtherPara.MaxValue,OtherPara.binsize(Id_index)+1);
       % [N,edgesNew,bin] = histcounts(YUse(1,:),edges); 
%          assignin('base', 'YUse', YUse(1,:));     
%          assignin('base', 'edges', edges);
% dd
        IniDist{Id_index}=N/sum(N);
%         
%                  assignin('base', 'IniDist', IniDist{Id_index});
% dd
       % Prefactor=[];
        

        for jjj=1:size(seqsUse,1)
            
       % IniDistTraj=zeros(1,size(N,2));
        
       % IniDistTraj(bin(jjj)+1)=IniDist{Id_index}(bin(jjj));
        
        %disp(IniDistTraj);disp(IniDist{Id_index});dd
        MatrixTemp2=cell(1,Nq);MatrixTemp2_2=cell(1,Nq);
        TransitionKernal_PartI=cell(1,Nq);TransitionKernal=cell(1,Nq);
        TrajEntropyFilteredTemp=cell(1,Nq);TransitionKernal_Product=cell(1,Nq);
        
          for Id_index_2=Id_index   
        MatrixTemp2{Id_index_2}=diag(ones(1,size(estTR{1,Id_index},1)));%/size(estTR{1,Id_index},1);
         %MatrixTemp2_2{Id_index_2}=diag(ones(1,size(estTR{1,Id_index},1)));
          end
          
          
        for kk=1:size(seqsUse,2)        
            for Id_index_2=Id_index
            
            %D = diag(estE{1,Id_index}(seqsUse(jjj,kk),:)); %previous 
            D = diag(estE{1,Id_index_2}(:,seqsUse(jjj,kk))); 
            %disp(size(D));disp(diag(D));ddd
            %MatrixTemp2=MatrixTemp2*estTR{1,Id_index}*D;
            %if kk>1
             MatrixTemp2{Id_index_2}=MatrixTemp2{Id_index_2}*estTR{1,Id_index_2}*D;
            %end
            % MatrixTemp2_2{Id_index_2}=MatrixTemp2_2{Id_index_2}*estTR{1,Id_index_2};% Point-wise
             %MatrixTemp2_2=MatrixTemp2_2*estTR{1,Id_index_2};
            %e = sort(real(eig(MatrixTemp2')),'descend');
            
            %[Q,R] = qr(MatrixTemp2');e=diag(R);%QR-decompoositon
            
            
% % % %             if Id_index_2==Id_index
% % % %             [Q_temp,R_temp] = qr(MatrixTemp2{Id_index_2}); R=R_temp'; %LQ-decompoositon
% % % %            
% % % %             e=diag(R);%e(e==0)=[];
% % % %             
% % % %             e=sort(abs(e),'descend'); %Keep zero or not, now not...
% % % %            %e=sort(e,'descend'); %Keep zero or not, now not...
% % % %            % e=sort(e,'descend'); %Keep zero or not, now not...
% % % % %             disp(e(1:5));
% % % % %             if kk>5
% % % % %                 dd
% % % % %             end
% % % %             if min(min(MatrixTemp2{Id_index_2}))<0
% % % %             disp(MatrixTemp2{Id_index_2});disp(e);dd
% % % %             end
% % % %             if e(1)<=0
% % % %                 count1=count1+1;
% % % %                 eTemp(1,kk)=realmin;
% % % %             else
% % % %                 eTemp(1,kk)=e(1);
% % % %             end
% % % %             if e(2)<=0
% % % %                 eTemp(2,kk)=realmin;count2=count2+1;
% % % %             else
% % % %                 eTemp(2,kk)=e(2);
% % % %             end
% % % %             %HMM.TrajEigen1{Id_index}(jjj,kk)=eTemp(1,kk);%e(1); 
% % % %             %HMM.TrajEigen2{Id_index}(jjj,kk)=eTemp(2,kk);%e(2);
% % % %             HMM.TrajLyapunov1{Id_index}(jjj,kk)=-log2(eTemp(1,kk))/kk; %sum(log2(eTemp(1,1:kk)+realmin))/kk; 
% % % %             HMM.TrajLyapunov2{Id_index}(jjj,kk)=-log2(eTemp(2,kk))/kk;%sum(log2(eTemp(2,1:kk)+realmin))/kk;
% % % %             %HMM.TrajLyapunov1{Id_index}(jjj,kk)=sum(log2(HMM.TrajEigen1{Id_index}(jjj,1:kk)+realmin))/kk; 
% % % %           %  HMM.TrajLyapunov2{Id_index}(jjj,kk)=sum(log2(HMM.TrajEigen2{Id_index}(jjj,1:kk)+realmin))/kk;
% % % %             end
            
%         if length(e)>2    
%         HMM.TrajEigen1{Id_index}(jjj,kk)=e(1); 
%         HMM.TrajEigen2{Id_index}(jjj,kk)=e(2);
%         HMM.TrajLyapunov1{Id_index}(jjj,kk)=sum(log(HMM.TrajEigen1{Id_index}(jjj,:)))/kk; 
%         HMM.TrajLyapunov2{Id_index}(jjj,kk)=sum(log(HMM.TrajEigen2{Id_index}(jjj,:)))/kk;
%         elseif length(e)>1
%             disp(R);disp(Id_index);dd
%         HMM.TrajEigen1{Id_index}(jjj,kk)=e(1);
%         HMM.TrajLyapunov1{Id_index}(jjj,kk)=sum(log(HMM.TrajEigen1{Id_index}(jjj,:)))/kk; 
%         HMM.TrajEigen2{Id_index}(jjj,kk)=realmin;
%         HMM.TrajLyapunov2{Id_index}(jjj,kk)=0;%HMM.TrajLyapunov2{Id_index}(jjj,kk-1);        
%         else
%          HMM.TrajEigen1{Id_index}(jjj,kk)=realmin; 
%          HMM.TrajLyapunov1{Id_index}(jjj,kk)=0;%HMM.TrajLyapunov1{Id_index}(jjj,kk-1);        
%         HMM.TrajEigen2{Id_index}(jjj,kk)=realmin;   
%          HMM.TrajLyapunov2{Id_index}(jjj,kk)=0;%HMM.TrajLyapunov2{Id_index}(jjj,kk-1);                
%         end
        
        
         %disp(size(IniDist{Id_index}));disp(size(MatrixTemp2))
        %MatrixTemp3=IniDist{Id_index}*pinv(estE{1,Id_index_2})*MatrixTemp2;%IniDist{Id_index}*MatrixTemp2;%IniDist{Id_index}*pinv(estE{1,Id_index_2})*MatrixTemp2;%dd
        %MatrixTemp3=IniDistTraj*pinv(estE{1,Id_index_2})*MatrixTemp2;%IniDist{Id_index}*MatrixTemp2;%IniDist{Id_index}*pinv(estE{1,Id_index_2})*MatrixTemp2;%dd
        
        if kk==1
            %IniDistHiddenState=estE{1,Id_index_2}(:,bin(jjj))'/sum(estE{1,Id_index_2}(:,bin(jjj)));
            %%Assumed initial hidden states distribution based on observed
            IniDistHiddenState=ones(1,size(estTR{1,Id_index},1))/sum(ones(1,size(estTR{1,Id_index},1)));
            %%Assumed uniform initial hidden states distribution 
            HMM.TrajProb{Id_index_2,Id_index}(jjj,1)=IniDist{Id_index}(bin(jjj));%IniDistTraj(bin(jjj)+1);%should be the observed probability           
           
                        %assignin('base', 'IniDistHiddenState', IniDistHiddenState);
            % assignin('base', 'MatrixTemp2', MatrixTemp2{Id_index_2});
%             assignin('base', 'IniDist', IniDist);
%             ddd
           % HMM.TrajProb_Point{Id_index_2,Id_index}(jjj,1)=IniDistTraj(bin(jjj)+1);%should be the observed probability
%             TrajEntropyFilteredTemp{1,Id_index_2}=-IniDistHiddenState.*log2(IniDistHiddenState+realmin);
%             HMM.TrajEntropyFiltered{Id_index_2,Id_index}(jjj,1)=sum(TrajEntropyFilteredTemp{1,Id_index_2});
            %HMM.TrajEntropyFiltered{Id_index_2,Id_index}(jjj,1)=abs(sum(IniDistHiddenState.*log2(IniDistHiddenState+realmin)));
            
              %assignin('base', 'IniDistHiddenState', IniDistHiddenState);ddd
%              assignin('base', 'TrajEntropyFilteredTemp', TrajEntropyFilteredTemp);
        end
        
        ProbabilityTemp3=IniDistHiddenState*MatrixTemp2{Id_index_2};
       % ProbabilityTemp3_2=IniDistHiddenState*MatrixTemp2_2{Id_index_2}*D;%new assumption
        HMM.TrajProb{Id_index_2,Id_index}(jjj,kk)=sum(ProbabilityTemp3);%sum(MatrixTemp3);%real(sum(MatrixTemp3));
        %HMM.TrajProb_Point{Id_index_2,Id_index}(jjj,kk+1)=sum(ProbabilityTemp3_2);%sum(MatrixTemp3);%real(sum(MatrixTemp3));
        %HMM.TrajProb{Id_index_2,Id_index}(jjj,kk+1)=sum(ProbabilityTemp3);%sum(MatrixTemp3);%real(sum(MatrixTemp3));
        %HMM.TrajProb_Point{Id_index_2,Id_index}(jjj,kk+1)=sum(ProbabilityTemp3_2);%sum(MatrixTemp3);%real(sum(MatrixTemp3));
        
%         
%         if sum(ProbabilityTemp3)>0
%             ProbabilityTemp3=ProbabilityTemp3/sum(ProbabilityTemp3);%filtered probability
%         else
%             ProbabilityTemp3=zeros(1,size(ProbabilityTemp3,2));%disp(ProbabilityTemp3);
%         end
%       
       
%         if kk==1
%             TransitionKernal_Product{1,Id_index_2}=diag(ones(1,size(estTR{1,Id_index},1)));
%         else
%             TransitionKernal_Product{1,Id_index_2}=TransitionKernal_Product{1,Id_index_2}*TransitionKernal{1,Id_index_2}; %multiply by the previous TransitionKernal
%         end
%         TransitionKernal_PartI{1,Id_index_2}=HMM.TrajProb{Id_index_2,Id_index}(jjj,kk)/HMM.TrajProb{Id_index_2,Id_index}(jjj,kk+1);
%         TransitionKernal{1,Id_index_2}=TransitionKernal_PartI{1,Id_index_2}*estTR{1,Id_index_2}*D;
%         
%         TrajEntropyFilteredTemp{1,Id_index_2}=TrajEntropyFilteredTemp{1,Id_index_2}*TransitionKernal{1,Id_index_2}...
%             -IniDistHiddenState*TransitionKernal_Product{1,Id_index_2}*(TransitionKernal{1,Id_index_2}.*log2(TransitionKernal{1,Id_index_2}+realmin));
%         HMM.TrajEntropyFiltered{Id_index_2,Id_index}(jjj,kk+1)=sum(TrajEntropyFilteredTemp{1,Id_index_2});
        
%          assignin('base', 'TransitionKernal_PartI', TrajEntropyFilteredTemp*TransitionKernal);
%          assignin('base', 'TransitionKernal', TransitionKernal);       
%          assignin('base', 'TrajEntropyFilteredTemp2', TrajEntropyFilteredTemp);
%          ddd
        %HMM.TrajEntropyFiltered{Id_index_2,Id_index}(jjj,kk+1)=abs(sum(ProbabilityTemp3.*log2(ProbabilityTemp3+realmin)));

        
        
        end
       % HMM.TrajMemory{Id_index}(jjj,kk)=abs(HMM.TrajLyapunov2{Id_index}(jjj,kk)-HMM.TrajLyapunov1{Id_index}(jjj,kk));
       
        end
        end
        
        %Remove Traj With NaN and zero
        kkk=0;
        for Id_index_2=Id_index
            DeleteRowsNaN(:,Id_index_2)=any(isnan(HMM.TrajProb{Id_index_2,Id_index}),2);
            DeleteRowsZero(:,Id_index_2)=any(HMM.TrajProb{Id_index_2,Id_index}==0,2);
        %HMM.TrajProb{Id_index_2,Id_index}(any(isnan(HMM.TrajProb{Id_index_2,Id_index}),2),:) = realmin;
        %HMM.TrajProb_Point{Id_index_2,Id_index}(any(isnan(HMM.TrajProb_Point{Id_index_2,Id_index}),2),:) = [];
%         HMM.TrajEntropyFiltered{Id_index_2,Id_index}(any(isnan(HMM.TrajEntropyFiltered{Id_index_2,Id_index}),2),:) = [];
            if kkk==0
                DeleteRowsNaNTotal=DeleteRowsNaN(:,Id_index_2);
                DeleteRowsZeroTotal=DeleteRowsZero(:,Id_index_2);
            else
                DeleteRowsNaNTotal=or(DeleteRowsNaNTotal,DeleteRowsNaN(:,Id_index_2));
                DeleteRowsZeroTotal=or(DeleteRowsZeroTotal,DeleteRowsZero(:,Id_index_2));
            end
            kkk=kkk+1;
        end
            DeleteRowsTotal=or(DeleteRowsNaNTotal,DeleteRowsZeroTotal);
            
        for Id_index_2=Id_index
        HMM.TrajProb{Id_index_2,Id_index}(DeleteRowsTotal,:) = [];%realmin;%[];
        %HMM.TrajProb{Id_index_2,Id_index}(any(HMM.TrajProb{Id_index_2,Id_index}==0,2),:) = [];%1e-10;realmin;%realmin;%[];
        end
        
       % HMM.TrajMemory{Id_index}(any(isnan(HMM.TrajMemory{Id_index}),2),:) = [];
        
        
%         if Id_index==1
%         OtherPara.CompensationFactor=0;
%         for Id_index_2=Id_index
%         OtherPara.CompensationFactor=OtherPara.CompensationFactor+size(Y{Id_index_2,1},2);%size(HMM.TrajProb{Id_index_2,1},1);
%         disp(OtherPara.CompensationFactor);
%         end       
%         OtherPara.CompensationFactor=OtherPara.CompensationFactor/Nq;
%         disp('CompensationFactor in HMM statistics');
%         disp(OtherPara.CompensationFactor);
%         end
        OtherPara.CompensationFactor=500; 
        
        for Id_index_2=Id_index
%         if OtherPara.TrajBasedNormalize==1
%          %   assignin('base', 'TrajProb1', HMM.TrajProb);
%         HMM.TrajProb{Id_index_2,Id_index}(:,:)=HMM.TrajProb{Id_index_2,Id_index}(:,:)./sum(HMM.TrajProb{Id_index_2,Id_index}(:,:),1);
% %           
%             
%             
%        % HMM.TrajProb_Point{Id_index_2,Id_index}(:,:)=HMM.TrajProb_Point{Id_index_2,Id_index}(:,:)./sum(HMM.TrajProb_Point{Id_index_2,Id_index}(:,:),1);
%        % HMM.TrajEntropy{Id_index_2,Id_index}=- (1/(kk+1)).*log2(HMM.TrajProb{Id_index_2,Id_index})/size(HMM.TrajProb{Id_index_2,Id_index},1)*OtherPara.CompensationFactor;
%        % HMM.TrajEntropy_Point{Id_index_2,Id_index}=- HMM.TrajProb_Point{Id_index_2,Id_index}.*log2(HMM.TrajProb_Point{Id_index_2,Id_index}+realmin)/size(HMM.TrajProb_Point{Id_index_2,Id_index},1)*OtherPara.CompensationFactor;
%         else 
%         
%          Prefactor=[];   
%          for nn=1:size(HMM.TrajProb{Id_index_2,Id_index},1)
%             Prefactor(nn,:)=1./[1:1:size(HMM.TrajProb{Id_index_2,Id_index},2)];
%          end
%         
%         HMM.TrajEntropy{Id_index_2,Id_index}=- Prefactor.*log2(HMM.TrajProb{Id_index_2,Id_index});
%         
%          %assignin('base', 'Prefactor', Prefactor);dddddd
%        % HMM.TrajEntropy_Point{Id_index_2,Id_index}=- HMM.TrajProb_Point{Id_index_2,Id_index}.*log2(HMM.TrajProb_Point{Id_index_2,Id_index}+realmin);
%         end
        
        
%         HMM.TrajEntropy{Id_index_2,Id_index}(any(isnan(HMM.TrajEntropy{Id_index_2,Id_index}),2),:) = [];
      %   HMM.TrajEntropy{Id_index_2,Id_index}(any(HMM.TrajEntropy{Id_index_2,Id_index}==Inf,2),:) = [];
         
%          %To save space
%          if Id_index_2~=Id_index
%              HMM.TrajEntropy{Id_index_2,Id_index}=[];
%          end
         
        end
        %HMM.TrajEntropy{Id_index}(any(isnan(HMM.TrajEntropy{Id_index}),2),:) = []; 
        
    
       
    end
    
    
%disp(count1);disp(count2);
DistanceConditions=0;
return;

    %%
    %plot ordered transition and emssion matrix
    figure ('position', [00, 10, 400, 800])
    
    for i=1:size(estE{1,Id_index},1)
        estSum(i)=sum(estE{1,Id_index}(i,:).*[1:size(estE{1,Id_index},2)]);
    end
    [B,estIndex] = sort(estSum);
    estEsort{1,Id_index}=estE{1,Id_index}(estIndex,:);
    estTRsort{1,Id_index}=estTR{1,Id_index}(estIndex,estIndex);
    %dd
%     
%     subplot(2,2,1)
%     
%     surf(estTR{1,Id_index}); title('Transition');set(gca,'FontSize',16); %shading interp;
%     view(2)
%    subplot(2,2,2)
%   
%    surf(estE{1,Id_index});  title('Emission'); set(gca,'FontSize',16);% shading interp;
%    view(2)
   subplot(2,1,1)
   
    surf(estTRsort{1,Id_index}); title('Ordered transition'); set(gca,'FontSize',16);%shading interp; 
    view(2)
   subplot(2,1,2)
  
   surf(estEsort{1,Id_index});  title('Ordered emission'); set(gca,'FontSize',16);%shading interp; 
   view(2)

    figurenamehmm=[foldername,'\Matrix_FitToAll_',num2str(ID(1)),'States_',num2str(OtherPara.state(Id_index)),'bin_',num2str(OtherPara.binsize(Id_index)),'.jpg'];
 saveas(gcf,figurenamehmm);

%end
  

 

end
 
 else
    
        YUse=seqSample_test'*OtherPara.Conversion(Id_index);
%         for jj=1:size(YUse,2)
%         [N,edgesNew,seqs2{Id_index}(jj,:)] = histcounts(YUse(:,jj)',OtherPara.edges{Id_index});
%          end

         seqsUse=seqSample_test;%seqs2{Id_index};
         [N,edgesNew,bin] = histcounts(YUse(1,:),OtherPara.edges{Id_index});
        IniDist{Id_index}=N/sum(N);
        
        for jjj=1:size(seqsUse,1)
 
        MatrixTemp2_2=cell(1,Nq);
          for Id_index_2=Id_index   
        MatrixTemp2_2{Id_index_2}=diag(ones(1,size(estTR{1,Id_index},1)));
          end
        for kk=1:size(seqsUse,2)        
            for Id_index_2=Id_index
            D_2 = diag(estE{1,Id_index_2}(:,seqsUse(jjj,kk))); 
             MatrixTemp2_2{Id_index_2}=MatrixTemp2_2{Id_index_2}*estTR{1,Id_index_2}*D_2;
   
        if kk==1
            IniDistHiddenState=ones(1,size(estTR{1,Id_index},1))/sum(ones(1,size(estTR{1,Id_index},1)));
            HMM.TrajProb2{Id_index_2,Id_index}(jjj,1)=IniDist{Id_index}(bin(jjj));%IniDistTraj(bin(jjj)+1);%should be the observed probability           
         end
        
        ProbabilityTemp3_2=IniDistHiddenState*MatrixTemp2_2{Id_index_2};
       HMM.TrajProb2{Id_index_2,Id_index}(jjj,kk)=sum(ProbabilityTemp3_2);%sum(MatrixTemp3);%real(sum(MatrixTemp3));      
            end
        
        end
        end
%           assignin('base', 'HMM', HMM);
%              assignin('base', 'seqSample_test', seqSample_test);
%               assignin('base', 'seqsUse', seqsUse);
%               assignin('base', 'seqs2', seqs2);
        kkk=0;
        for Id_index_2=Id_index
            DeleteRowsNaN(:,Id_index_2)=any(isnan(HMM.TrajProb2{Id_index_2,Id_index}),2);
            DeleteRowsZero(:,Id_index_2)=any(HMM.TrajProb2{Id_index_2,Id_index}==0,2);
            if kkk==0
                DeleteRowsNaNTotal=DeleteRowsNaN(:,Id_index_2);
                DeleteRowsZeroTotal=DeleteRowsZero(:,Id_index_2);
            else
                DeleteRowsNaNTotal=or(DeleteRowsNaNTotal,DeleteRowsNaN(:,Id_index_2));
                DeleteRowsZeroTotal=or(DeleteRowsZeroTotal,DeleteRowsZero(:,Id_index_2));
            end
            kkk=kkk+1;
        end
            DeleteRowsTotal=or(DeleteRowsNaNTotal,DeleteRowsZeroTotal);
        
          
        for Id_index_2=Id_index
        HMM.TrajProb2{Id_index_2,Id_index}(DeleteRowsTotal,:) = [];%realmin;%[];
       end
       
        
        
end
 end
 close all;

end
disp('Calculate Statistics done!');

%  if OtherPara.HMMFitToAll==1
% disp('Calculate Channel capacity');     
% %Get the optimal trajectory mutual information
% %Need do after the whole loop to get statistics for all conditions
% options = optimoptions('fmincon','Algorithm','active-set');%,'Display','iter');
% % Define information calculation function (getInfo.m), initialize bounds for optimization.
% %infoF= @(q) -getInfo(q,fX1_all,fX_all,timeunit);
% %MI_Full= @(q) -getInfoTrajHMM(q,HMM);
% Nq1=Nq;
% %HMM = rmfield(HMM,'MI_Full');
% 
% 
% if OtherPara.PointWiseOptim==1%OtherPara.PointWiseOptim==1
%     %return;
%     %assignin('base', 'MI_Full', MI_Full);
%     for m=1:size(HMM.TrajProb{1,1},2)
%         disp(m);
%         %MI_FullTarget=MI_Full(m);
%         OtherPara.Point=m;
%         MI_Full1= @(q) -getInfoHMM_New(q,HMM,Nq1,OtherPara);
%         Aeq=ones(1,Nq1);
%         beq=1;
%         LB=ones(Nq1,1)*realmin;
%         UB=ones(Nq1,1);  
%         Qstart=ones(Nq1,1);
%         Qstart=Qstart/sum(Qstart); 
%         HMM.qfit=fmincon(MI_Full1,Qstart,[],[],Aeq,beq,LB,UB,[],options);
%         %disp(HMM.qfit);
%         %[info1,HMM.MI_Full]=getInfoTrajHMM(HMM.qfit,HMM);
%         [info1,MI_Full1]=getInfoHMM_New(HMM.qfit,HMM,Nq1,OtherPara);        
%         HMM.MI_Full(m)=MI_Full1;
% %         if m==17
% %             assignin('base', 'HMMMI_Full', HMM.MI_Full);
% %             ddd
% %         end
%     end
% 
%     
%     
% else
%    
%     
% MI_Full= @(q) -getInfoHMM_Matrice(q,HMM,Nq1,OtherPara);
% Aeq=ones(1,Nq1);
% beq=1;
% LB=ones(Nq1,1)*realmin;
% UB=ones(Nq1,1);    
% Qstart=ones(Nq1,1);
% Qstart=Qstart/sum(Qstart); 
% HMM.qfit=fmincon(MI_Full,Qstart,[],[],Aeq,beq,LB,UB,[],options);
% disp(HMM.qfit);
% %[info1,HMM.MI_Full]=getInfoTrajHMM(HMM.qfit,HMM);
% [info1,HMM.MI_Full]=getInfoHMM_Matrice(HMM.qfit,HMM,Nq1,OtherPara);
% end
% 
% 
% 
%  end