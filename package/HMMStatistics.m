% Calculate trajectory probability for each time series from the hidden
% Markov model. Then, calculate the trajectory entropy for each condition, and time-dependent
% channel capacity for the dataset.

function [DistanceConditions,HMM]=HMMStatistics(estTR,estE,OtherPara,seqs,Y,Nq,foldername,ID,HMM,SampledTraj)

%OtherPara.Jackknife=0;
DistanceConditions=0;
    TrajSample=cell(1,Nq);
    TrajSample_Smooth=cell(1,Nq);
    TrajSample_HiddenStates=cell(1,Nq);
    estTRsort=cell(1,Nq);
    estEsort=cell(1,Nq);
    DistanceTrajs=cell(1,Nq);Distance2_Trajs=cell(1,Nq);
    
    
    if OtherPara.TrajBasedStatistic==1 
    HMM.TrajProb=cell(Nq,Nq);HMM.TrajProb_Point=cell(Nq,Nq);
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
    
    %Calculate model/condition based statistics
if OtherPara.HMMFitToAll==1 && OtherPara.CalculateStatistics==1 && OtherPara.TrajBasedStatistic==1 
    %disp(Id_index);
    HMM.on=0;
    if OtherPara.ModelBasedStatistic==1
        [V,D] = eig(estTR{1,Id_index}');
        e= real(eig(estTR{1,Id_index}'));%-eye(size(estTR{1,Id_index},1)));
        [M,I]=min(abs(e-1));
        X=abs(real(V(:,I)));X=X/sum(X);

        MatrixTemp1=estTR{1,Id_index};

        HMM.Entropy(Id_index)=0;HMM.Entropy_2(Id_index)=0;
        HMM.EntropyProduction(Id_index)=0;HMM.EntropyProduction_2(Id_index)=0;HMM.EntropyProduction00(Id_index)=0;

         EntropyProductionTest=NaN(size(MatrixTemp1,1),size(MatrixTemp1,2));
        for ii=1:size(MatrixTemp1,1)
        for jj=1:size(MatrixTemp1,2)
            if MatrixTemp1(ii,jj)>0
                HMM.Entropy(Id_index)=HMM.Entropy(Id_index)-X(ii)*MatrixTemp1(ii,jj)*log2(MatrixTemp1(ii,jj)+realmin);
                if MatrixTemp1(jj,ii)>0
                    
                EntropyProductionTest(ii,jj)=MatrixTemp1(ii,jj)/MatrixTemp1(jj,ii);
                if X(ii)*MatrixTemp1(ii,jj)*log2(X(ii)*MatrixTemp1(ii,jj)/(X(jj)*MatrixTemp1(jj,ii)))<Inf
                HMM.EntropyProduction(Id_index)=HMM.EntropyProduction(Id_index)+X(ii)*MatrixTemp1(ii,jj)...
                    *log2(X(ii)*MatrixTemp1(ii,jj)/(X(jj)*MatrixTemp1(jj,ii)));
                end          
                HMM.EntropyProduction00(Id_index)=HMM.EntropyProduction00(Id_index)+X(ii)*MatrixTemp1(ii,jj)...
                    *log2(MatrixTemp1(ii,jj)/(MatrixTemp1(jj,ii)+realmin)+realmin);
                end
            end
        end
        end

    end
    
  %  continue;
    %%
     %Calculate traj based statistics
    if OtherPara.TrajBasedStatistic==1
        
         DeleteRowsNaN=[];DeleteRowsZero=[];
            
            
        if Id_index==Nq+1
         seqsUse=cat(1,seqs{1:Nq});
         YUse=cat(2,Y{1:Nq});
        else
          seqsUse=seqs{Id_index};
         YUse=Y{Id_index};
        end
        
        edges=linspace(OtherPara.MinValue,OtherPara.MaxValue,OtherPara.binsize+1);
        [N,edgesNew,bin] = histcounts(YUse(1,:),edges); 

        IniDist{Id_index}=N/sum(N);
%         


        for jjj=1:min(size(seqsUse,1),size(YUse,1))
            
        MatrixTemp2=cell(1,Nq);MatrixTemp2_2=cell(1,Nq);
%         TransitionKernal_PartI=cell(1,Nq);TransitionKernal=cell(1,Nq);
%         TrajEntropyFilteredTemp=cell(1,Nq);TransitionKernal_Product=cell(1,Nq);
        
          for Id_index_2=1:Nq   
        MatrixTemp2{Id_index_2}=diag(ones(1,size(estTR{1,Id_index},1)));%/size(estTR{1,Id_index},1);
         %MatrixTemp2_2{Id_index_2}=diag(ones(1,size(estTR{1,Id_index},1)));
          end
          
          
        for kk=1:size(seqsUse,2)        
            for Id_index_2=1:Nq
            
            %D = diag(estE{1,Id_index}(seqsUse(jjj,kk),:)); %previous 
            D = diag(estE{1,Id_index_2}(:,seqsUse(jjj,kk))); 
             MatrixTemp2{Id_index_2}=MatrixTemp2{Id_index_2}*estTR{1,Id_index_2}*D;
         
            
            
            if Id_index_2==Id_index
            [Q_temp,R_temp] = qr(MatrixTemp2{Id_index_2}); R=R_temp'; %LQ-decompoositon
           
            e=diag(R);%e(e==0)=[];
            
            e=sort(abs(e),'descend'); 
            if min(min(MatrixTemp2{Id_index_2}))<0
            disp(MatrixTemp2{Id_index_2});disp(e);dd
            end
            if e(1)<=0
                count1=count1+1;
                eTemp(1,kk)=realmin;
            else
                eTemp(1,kk)=e(1);
            end
            if e(2)<=0
                eTemp(2,kk)=realmin;count2=count2+1;
            else
                eTemp(2,kk)=e(2);
            end
            HMM.TrajLyapunov1{Id_index}(jjj,kk)=-log2(eTemp(1,kk))/kk; %sum(log2(eTemp(1,1:kk)+realmin))/kk; 
            HMM.TrajLyapunov2{Id_index}(jjj,kk)=-log2(eTemp(2,kk))/kk;%sum(log2(eTemp(2,1:kk)+realmin))
            end
            
       
        if kk==1
            IniDistHiddenState=ones(1,size(estE{1,Id_index_2},1))/size(estE{1,Id_index_2},1);%estE{1,Id_index_2}(:,bin(jjj))'/sum(estE{1,Id_index_2}(:,bin(jjj)));%new assumption
            %disp(Id_index);disp(Id_index_2);disp(jjj);disp(size(bin));disp(size(IniDist{Id_index}));
            HMM.TrajProb{Id_index_2,Id_index}(jjj,1)=IniDist{Id_index}(bin(jjj));%IniDistTraj(bin(jjj)+1);%should be the observed probability           
            
        end
        
        ProbabilityTemp3=IniDistHiddenState*MatrixTemp2{Id_index_2};
        HMM.TrajProb{Id_index_2,Id_index}(jjj,kk)=sum(ProbabilityTemp3);%sum(MatrixTemp3);%real(sum(MatrixTemp3));
    
        end
        HMM.TrajMemory{Id_index}(jjj,kk)=abs(HMM.TrajLyapunov2{Id_index}(jjj,kk)-HMM.TrajLyapunov1{Id_index}(jjj,kk));
       
        end
        end
      %   assignin('base', 'IniDistHiddenState', IniDistHiddenState);
      %  assignin('base', 'TrajProb1', HMM.TrajProb);
        
        %Remove Traj With NaN and zero
        for Id_index_2=1:Nq
            DeleteRowsNaN(:,Id_index_2)=any(isnan(HMM.TrajProb{Id_index_2,Id_index}),2);
            DeleteRowsZero(:,Id_index_2)=any(HMM.TrajProb{Id_index_2,Id_index}==0,2);
       
            if Id_index_2==1
                DeleteRowsNaNTotal=DeleteRowsNaN(:,Id_index_2);
                DeleteRowsZeroTotal=DeleteRowsZero(:,Id_index_2);
            else
                DeleteRowsNaNTotal=or(DeleteRowsNaNTotal,DeleteRowsNaN(:,Id_index_2));
                DeleteRowsZeroTotal=or(DeleteRowsZeroTotal,DeleteRowsZero(:,Id_index_2));
            end
        end
            DeleteRowsTotal=or(DeleteRowsNaNTotal,DeleteRowsZeroTotal);
            
        for Id_index_2=1:Nq
        HMM.TrajProb{Id_index_2,Id_index}(DeleteRowsTotal,:) = [];%realmin;%[];
        %HMM.TrajProb{Id_index_2,Id_index}(any(HMM.TrajProb{Id_index_2,Id_index}==0,2),:) = [];%1e-10;realmin;%realmin;%[];
        end
       % assignin('base', 'DeleteRowsTotal', DeleteRowsTotal);
        HMM.TrajMemory{Id_index}(any(isnan(HMM.TrajMemory{Id_index}),2),:) = [];
        %assignin('base', 'TrajProb2', HMM.TrajProb);

        OtherPara.CompensationFactor=500; 
        
%         for Id_index_2=1:Nq
       % if OtherPara.TrajBasedNormalize==1
         %   assignin('base', 'TrajProb1', HMM.TrajProb);
        %HMM.TrajProb{Id_index_2,Id_index}(:,:)=HMM.TrajProb{Id_index_2,Id_index}(:,:)./sum(HMM.TrajProb{Id_index_2,Id_index}(:,:),1);
        %else 
        
%          Prefactor=[];   
%          for nn=1:size(HMM.TrajProb{Id_index_2,Id_index},1)
%              
%             Prefactor(nn,:)=1./[1:1:size(HMM.TrajProb{Id_index_2,Id_index},2)];
%          end
%         
%         HMM.TrajEntropy{Id_index_2,Id_index}=- Prefactor.*log2(HMM.TrajProb{Id_index_2,Id_index});
%         
%         %end
%         
%         
%          HMM.TrajEntropy{Id_index_2,Id_index}(any(isnan(HMM.TrajEntropy{Id_index_2,Id_index}),2),:) = [];
%          HMM.TrajEntropy{Id_index_2,Id_index}(any(HMM.TrajEntropy{Id_index_2,Id_index}==Inf,2),:) = [];
%          
%          %To save space
%          if Id_index_2~=Id_index
%              HMM.TrajEntropy{Id_index_2,Id_index}=[];
%          end
%          
%         end
        
    
       
    end
    

DistanceConditions=0;
%continue;

    %%
    %plot ordered transition and emssion matrix
     %('position', [00, 10, 400, 800])
    
    for i=1:size(estE{1,Id_index},1)
        estSum(i)=sum(estE{1,Id_index}(i,:).*[1:size(estE{1,Id_index},2)]);
    end
    [B,estIndex] = sort(estSum);
    estEsort{1,Id_index}=estE{1,Id_index}(estIndex,:);
    estTRsort{1,Id_index}=estTR{1,Id_index}(estIndex,estIndex);
%    subplot(2,1,1)   
%     surf(estTRsort{1,Id_index}); title('Ordered transition'); set(gca,'FontSize',16);%shading interp; 
%     view(2)
    myFig=figure('position', [00, 10, 600, 600]);
    surfc(estTRsort{1,Id_index});
    h = colorbar;
    ylabel(h, 'Probability')
    caxis([0 1]);
    colormap jet;
    view(2)
    title('Transition matrix');
    xlabel('Post hidden state');
     set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    %title(['Time: ',num2str(round(k/143*12,2)),' h'])
    ylabel('Pre hidden state')
     set(gca,'FontSize',40);
    set(findall(myFig, 'Type', 'Text'),'FontWeight', 'Normal')
    figurenamehmm=[foldername,'\TransitionMatrix_',num2str(ID(Id_index)),'States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize),'.jpg'];
 saveas(gcf,figurenamehmm);
%    subplot(2,1,2)
%    surf(estEsort{1,Id_index});  title('Ordered emission'); set(gca,'FontSize',16);%shading interp; 
%    view(2)

myFig=figure('position', [00, 10, 600, 600]);
 surfc(estEsort{1,Id_index});
    h = colorbar;
    ylabel(h, 'Probability')
    caxis([0 1]);
    colormap jet;
    view(2)
    title('Emission matrix');
    xlabel('Emission state');
     set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    %title(['Time: ',num2str(round(k/143*12,2)),' h'])
    ylabel('Hidden state')
    set(gca,'FontSize',40);
    set(findall(myFig, 'Type', 'Text'),'FontWeight', 'Normal')
    
    figurenamehmm=[foldername,'\EmissionMatrix_',num2str(ID(Id_index)),'States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize),'.jpg'];
 saveas(gcf,figurenamehmm);

%end
  

 

    end
 
 
    close all;

end
disp('Calculate Statistics done!');

 if OtherPara.HMMFitToAll==1
disp('Calculate Channel capacity');     
%Get the optimal trajectory mutual information
%Need do after the whole loop to get statistics for all conditions

%options = optimoptions('fmincon','Algorithm','active-set');%,'Display','iter');
options = optimoptions('fmincon','Algorithm','interior-point','Display','off');
%options = optimoptions('fmincon','Algorithm','sqp');%,'Display','iter');


Nq1=Nq;
if OtherPara.CCSelect==1 %&& Nq1>10%calculate CC for selected conditions, such as pairwise or 4 conditions etc   
    Nq1=length(OtherPara.CCSelectCondi);
end

%HMM = rmfield(HMM,'MI_Full');

cutoff=realmin;

   

if OtherPara.CCSelect~=1% For selected conditions, do not do Jackknife
    HMM.tempTrajProb=cell(Nq,Nq);ktemp=0;
    for kkk=OtherPara.JackknifeRatio
        ktemp=ktemp+1;
        disp('Jackknife procedure');

        for ll=1:size(HMM.TrajProb,2)
            remain=round(kkk*size(HMM.TrajProb{1,ll},1));
            for mm=1:size(HMM.TrajProb,1)
                %HMM.TrajProb{mm,ll}(remain:end,:)=[];
                %HMM.tempTrajProb{mm,ll}=HMM.TrajProb{mm,ll}(1:remain,:);%HMM.TrajProb{mm,ll};
                HMM.tempTrajProb{mm,ll}=HMM.TrajProb{mm,ll}(end-remain+1:end,:);%HMM.TrajProb{mm,ll};
            end

        end
        if OtherPara.PointWiseOptim==1%1%OtherPara.PointWiseOptim==1
            disp('Pointwise optimization');
            [MI_Full qfit]=PointWiseOptim(HMM,Nq1,OtherPara,options);    
            HMM.qfit1(:,ktemp)=qfit; 
            
%             assignin('base', 'aaaa', qfit);
%              assignin('base', 'HMM2', HMM);ddd
            
            HMM.MI_Jackknife(ktemp,:)=MI_Full;
        else
            OtherPara.Point=size(HMM.tempTrajProb{1,1},2);
            %MI_Full1= @(q) -getInfoHMM_New(q,HMM,Nq1,OtherPara);
            MI_Full1= @(q) -getInfoHMM_NewSum(q,HMM,Nq1,OtherPara);
            Aeq=ones(1,Nq1);
            beq=1;
            LB=ones(Nq1,1)*cutoff;
            UB=ones(Nq1,1)-cutoff;
            Qstart=ones(Nq1,1);
            Qstart=Qstart/sum(Qstart);
            qfit=fmincon(MI_Full1,Qstart,[],[],Aeq,beq,LB,UB,[],options);
            HMM.qfit(:,ktemp)=qfit;
            [MI_Full1 MI_Fulltemp]=getInfoHMM_NewSum(HMM.qfit(:,ktemp),HMM,Nq1,OtherPara);
            %assignin('base', 'kk', MI_Fulltemp);
            HMM.MI_Jackknife(ktemp,:)=MI_Fulltemp;
        end
        

        %assignin('base', 'aaaa', yfit);
    end

        for kkk=1:size(HMM.MI_Jackknife,1)
            x(kkk)=1/OtherPara.JackknifeRatio(kkk);%y1(kk)=HMM.MI_Full(kk,end);
            y1(kkk)=mean(HMM.MI_Jackknife(kkk,:),2);
        end
        P = polyfit(x,y1,1);
        xx=0:0.1:x(end);
        yfit = P(1)*xx+P(2);
        %plot(xx,yfit,'--','color',[0.5 .5 .5],'linewidth',3);
        HMM.yIntersection=yfit(1);
        HMM.IntersectionRatio=HMM.yIntersection/y1(1);
end
    
    OtherPara.Jackknife=0;
    if OtherPara.PointWiseOptim==1%1%OtherPara.PointWiseOptim==1
        disp('Pointwise optimization');
        [MI_Full qfit]=PointWiseOptim(HMM,Nq1,OtherPara,options);    
        HMM.qfit=qfit;        
        HMM.MI_Full=max(MI_Full*HMM.IntersectionRatio,0);
        
       
    else
        disp('No-pointwise optimization');
        cutoff=realmin;
        OtherPara.Point=size(HMM.TrajProb{1,1},2);
        %MI_Full1= @(q) -getInfoHMM_New(q,HMM,Nq1,OtherPara);
        MI_Full1= @(q) -getInfoHMM_NewSum(q,HMM,Nq1,OtherPara);
        Aeq=ones(1,Nq1);
        beq=1;
        LB=ones(Nq1,1)*cutoff;
        UB=ones(Nq1,1)-cutoff;  
        Qstart=ones(Nq1,1);
        Qstart=Qstart/sum(Qstart); 
        HMM.qfit=fmincon(MI_Full1,Qstart,[],[],Aeq,beq,LB,UB,[],options);
        disp(HMM.qfit);
        % for m=1:size(HMM.TrajProb{1,1},2)
        %     OtherPara.Point=m;
        % MI_Full1=getInfoHMM_New(HMM.qfit,HMM,Nq1,OtherPara);        
        % HMM.MI_Full(m)=MI_Full1;
        % end
        [MI_Full1 HMM.MI_Full2]=getInfoHMM_NewSum(HMM.qfit,HMM,Nq1,OtherPara); 
        %[MI_Full qfit]=PointWiseOptim(HMM,Nq1,OtherPara,options);           
        HMM.MI_Full=MI_Full1*HMM.IntersectionRatio;
    end

    





 end