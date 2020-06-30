function [HMM]=MMStatistics(fX_all,OtherPara,Y,seqs,Nq)

%OtherPara.Jackknife=0;

   
   HMM.TrajProb=cell(Nq,Nq);
  


IniDist=cell(Nq,1);
for Id_index=1:Nq%+1
%disp(Id_index);
    
    HMM.on=0;
    
     %Calculate traj based statistics
   % if OtherPara.TrajBasedStatistic==1
        
         DeleteRowsNaN=[];DeleteRowsZero=[];
            
            
       
          seqsUse=seqs{Id_index};%seqs{Id_index};
         YUse=Y{Id_index};
   
        
        % assignin('base', 'YUse', YUse); dddd
        edges=linspace(OtherPara.MinValue,OtherPara.MaxValue,OtherPara.binsize+1);
        [N,edgesNew,bin] = histcounts(YUse(1,:),edges); 

        IniDist{Id_index}=N/sum(N);
%         


        for jjj=1:size(seqsUse,1) %for each traj
            
        %MatrixTemp2=cell(1,Nq);
        
          for Id_index_2=1:Nq   
        %MatrixTemp2{Id_index_2}=diag(ones(1,size(estTR{1,Id_index},1)));%/size(estTR{1,Id_index},1);
         %MatrixTemp2_2{Id_index_2}=diag(ones(1,size(estTR{1,Id_index},1)));
          end
          
          
        for kk=1:size(seqsUse,2)        
            for Id_index_2=1:Nq
            
            %D = diag(estE{1,Id_index}(seqsUse(jjj,kk),:)); %previous 
           % D = diag(estE{1,Id_index_2}(:,seqsUse(jjj,kk))); 
             %MatrixTemp2{Id_index_2}=MatrixTemp2{Id_index_2}*estTR{1,Id_index_2}*D;
               
       
            if kk==1
            %IniDistHiddenState=estE{1,Id_index_2}(:,bin(jjj))'/sum(estE{1,Id_index_2}(:,bin(jjj)));%new assumption
            HMM.TrajProb{Id_index_2,Id_index}(jjj,1)=IniDist{Id_index}(bin(jjj));%IniDistTraj(bin(jjj)+1);%should be the observed probability           

            else
        
            %ProbabilityTemp3=IniDistHiddenState*MatrixTemp2{Id_index_2};
            
            Joint=fX_all{Id_index_2,kk-1};
            Marginal=sum(Joint,2);
            transition_probabilities=Joint;
            transition_probabilities(Marginal>0,:)=Joint(Marginal>0,:)./Marginal(Marginal>0);
            
            HMM.TrajProb{Id_index_2,Id_index}(jjj,kk)=HMM.TrajProb{Id_index_2,Id_index}(jjj,kk-1)*transition_probabilities(seqsUse(jjj,kk-1),seqsUse(jjj,kk));%sum(MatrixTemp3);%real(sum(MatrixTemp3));
            
            %HMM.TrajProb{Id_index_2,Id_index}(jjj,kk)=HMM.TrajProb{Id_index_2,Id_index}(jjj,kk-1)*fX_all{Id_index_2,kk-1}(seqsUse(jjj,kk-1),seqsUse(jjj,kk));%sum(MatrixTemp3);%real(sum(MatrixTemp3));
    
            end
        %HMM.TrajMemory{Id_index}(jjj,kk)=abs(HMM.TrajLyapunov2{Id_index}(jjj,kk)-HMM.TrajLyapunov1{Id_index}(jjj,kk));
       
        end
        end
       

        end
        
        
         if OtherPara.PartialTrainTest~=1 %not do for likelihood calculation
        %Remove Traj With NaN and zero
        for Id_index_2=1:Nq
%             %DeleteRowsNaN(:,Id_index_2)=any(isnan(HMM.TrajProb{Id_index_2,Id_index}),2);
%             DeleteRowsZero(:,Id_index_2)=any(HMM.TrajProb{Id_index_2,Id_index}==0,2);
%        
%             if Id_index_2==1
%                 %DeleteRowsNaNTotal=DeleteRowsNaN(:,Id_index_2);
%                 DeleteRowsZeroTotal=DeleteRowsZero(:,Id_index_2);
%             else
%                 %DeleteRowsNaNTotal=or(DeleteRowsNaNTotal,DeleteRowsNaN(:,Id_index_2));
%                 DeleteRowsZeroTotal=or(DeleteRowsZeroTotal,DeleteRowsZero(:,Id_index_2));
%             end
            if Id_index_2==Id_index
                DeleteRowsZeroTotal=any(HMM.TrajProb{Id_index_2,Id_index}==0,2);
               % disp(find(DeleteRowsZeroTotal~=0))
            end
        end
        
        
           % DeleteRowsTotal=or(DeleteRowsNaNTotal,DeleteRowsZeroTotal);
            DeleteRowsTotal=DeleteRowsZeroTotal;
            
        for Id_index_2=1:Nq
        HMM.TrajProb{Id_index_2,Id_index}(DeleteRowsTotal,:) = [];%realmin;%[];
        %HMM.TrajProb{Id_index_2,Id_index}(any(HMM.TrajProb{Id_index_2,Id_index}==0,2),:) = [];%1e-10;realmin;%realmin;%[];
        end
         end
        
        
end


disp('Calculating Statistics done!');

if OtherPara.PartialTrainTest==1
    return;
end

%OtherPara.HMMFitToAll=1;OtherPara.Jackknife=1;OtherPara.PointWiseOptim=1;OtherPara.CCSelect=0;

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

       % assignin('base', 'aaaa', yfit);
    


%else
    
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