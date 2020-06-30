function [filenamehmm HMM]=VectorMethod(Y,OtherPara,MainPara)

LabelName=['VectorMethod_',OtherPara.LabelName];%'Brooks_HMM_Condition_All_2';
if OtherPara.HMM==4 %Time-point method
    LabelName=['TimePointMethod_',OtherPara.LabelName];%'Brooks_HMM_Condition_All_2';
end
        
% if MainPara.permuteData==1
% LabelName=[LabelName,'_Permute'];
% end

if MainPara.CCSelect==1
    if length(MainPara.CCSelectCondi)>length(Y) || max(MainPara.CCSelectCondi)>length(Y)
        disp('The seletive conditions index is out of range.');
        return;
    else
        YY=Y(MainPara.CCSelectCondi);
        Y=YY;
    end
end


filenamehmm=[pwd,'\',LabelName,'_cond',num2str(length(Y)),'.mat'];
if MainPara.permuteData==1
        filenamehmm=[pwd,'\',LabelName,'_cond',num2str(length(Y)),'_permute_',num2str(MainPara.permuteMode),'.mat'];    
end
if OtherPara.HMM~=4
if MainPara.WindowMode==2
    filenamehmm=[pwd,'\',LabelName,'_cond',num2str(length(Y)),'_Window_',num2str(MainPara.WindowLength),'.mat'];   
end
end

if exist(filenamehmm) && MainPara.CCSelect~=1
    disp('Have conducted the mutual information calculation.');
    disp(['filename ',filenamehmm]);
    load(filenamehmm);
    return;
end

tic;
%assignin('base', 'XX', Y);ddd
for TimePointsIndex=1:1:size(Y{1},1)%[1:5:size(txt,1)];
    disp(TimePointsIndex);toc;
    
     %Naive scheme
        if OtherPara.HMM~=4
        if MainPara.WindowMode==1
            %Naive scheme
           featureNames = 1:TimePointsIndex;% X{1,j}=data.nfkb(jj).metrics.time_series(:,1:TimePointsIndex)';
        elseif MainPara.WindowMode==2
            %Intermidate scheme: with using 5 time points
           
         Left=max(1,TimePointsIndex-floor(MainPara.WindowLength/2));
         Right=min(size(Y{1},1),TimePointsIndex+floor(MainPara.WindowLength/2));
         disp(Left);disp(Right);
         featureNames = Left:Right;%X{1,j}=data.nfkb(jj).metrics.time_series(:,Left:Right)';
        end
        end
        
        if OtherPara.HMM==4 %Time-point method
            featureNames=TimePointsIndex;
        end
        
        for k=1:length(Y)
            X{k}=Y{k}(featureNames,:);
        end

    if MainPara.permuteData==1
        [index X]=PermData(X,MainPara.permuteMode);        
    end
    
    
[fitI,info1,I,Q]= jacknifeCC(X,9,12); 
     %CC(TimePointsIndex)=info1;
     CC(TimePointsIndex)=max(info1,0);
end


HMM.MI_Full=CC;

save(filenamehmm,'HMM','OtherPara','MainPara');
end