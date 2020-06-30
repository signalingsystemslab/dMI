function [MI_Full qfit]=PointWiseOptim(HMM,Nq1,OtherPara,options)

% if OtherPara.CCSelect==1 && Nq1>20%calculate CC for selected conditions, such as pairwise or 4 conditions etc   
%     Nq1=length(OtherPara.CCSelectCondi);
% end



for m=1:size(HMM.TrajProb{1,1},2)
        if OtherPara.Jackknife==0
        disp(['Maximum mutual information at the time point ', num2str(m)]);
        end
        %MI_FullTarget=MI_Full(m);
        OtherPara.Point=m;
        %disp(Nq1);dddd
        MI_Full1= @(q) -getInfoHMM_New(q,HMM,Nq1,OtherPara);
        %MI_Full1= @(q) -getInfoHMM_NewSum(q,HMM,Nq1,OtherPara);
        Aeq=ones(1,Nq1);
        beq=1;
        LB=ones(Nq1,1)*0;%ones(Nq1,1)*eps;%ones(Nq1,1)*realmin;
        UB=ones(Nq1,1);  
        Qstart=ones(Nq1,1);
        Qstart=Qstart/sum(Qstart); 
        %nonlcon = @(x) -Constrain(x,HMM,Nq1,OtherPara);
        %qfit=fmincon(MI_Full1,Qstart,[],[],Aeq,beq,LB,UB,nonlcon,options);
        qfit=fmincon(MI_Full1,Qstart,[],[],Aeq,beq,LB,UB,[],options);
        %disp(HMM.qfit);
        %[info1,HMM.MI_Full]=getInfoTrajHMM(HMM.qfit,HMM);
        [MI_Full1Opt aSummaryTemp]=getInfoHMM_New(qfit,HMM,Nq1,OtherPara);        
        MI_Full(m)=max(MI_Full1Opt,0);
        
        aSummary.a0(:,m)=aSummaryTemp.a0(:,m);
        aSummary.a1(m)=aSummaryTemp.a1(m);
        aSummary.a2(:,m)=aSummaryTemp.a2(:,m);
        aSummary.a3(m)=aSummaryTemp.a3(m);
         aSummary.q(:,m)=qfit;%aSummaryTemp.q;
%         if m==95
%             assignin('base', 'HMMMI_Full', HMM.MI_Full);
%            % assignin('base', 'HMMMI_Full', HMM.MI_Full);
%             ddd
%         end
end
    %assignin('base', 'aSummary', aSummary);


end