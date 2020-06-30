function [index X]=PermData(X,permuteMode)

       % assignin('base', 'X', X);
       
       if permuteMode>3%concentrate cell
           XFull=[];
           for i=1:length(X)
               XFull=[XFull,X{i}];
           end 
           XFullVector=reshape(XFull,[],1);
           XFullVectorTemp = XFullVector(~isnan(XFullVector));
           XFullVector = XFullVectorTemp(isfinite(XFullVectorTemp));
           %isfinite(A)
           %disp(size(XFullVector));disp(size(XFullVectorTemp));
           index=0;
       end
       
        for i=1:length(X)
            
            if permuteMode==4%Brute force permutation
                disp('Random take points from all data.')
                for kk=1:size(X{i},2)
                Z{i}(:,kk)=datasample(XFullVector,size(X{i},1),'Replace',false);
                end 
                if i==2%length(X)
                    X=Z;break;
                end
                continue;
            end
            
            if permuteMode==5%Brute force permutation
                disp('Random take trajectories.')
                Z{i}=XFull(:,randi([1 size(XFull,2)],1,size(X{i},2))); 
                if i==2%length(X)
                    X=Z;break;
                end
                continue;
            end
            
            
            if permuteMode==6%Brute force permutation
                disp('Random take trajectories with adding fixed noise.')
                if i==1
                    Z{i}=max(XFull(:,randi([1 size(XFull,2)],1,size(X{i},2))),0);  
                elseif i==2
                    TempVector=reshape(X{1},[],1);
                    TempVector(isnan(TempVector))=[];
                    %Y{i}=max(Y{1}+std(TempVector)*0.02*randn(size(Y{1},1),size(Y{1},2)),0);
                    Z{i}=max(Z{1}+std(TempVector)*0.1*randn(size(Z{1},1),size(Z{1},2)),eps);
                    disp(std(TempVector));
                    X=Z;
                    break;
                end
                continue;
            end
            
            
            if permuteMode==7%Brute force permutation
                disp('Random take points from all data, same number of trajs.')
                for kk=1:size(X{1},2)
                Z{i}(:,kk)=datasample(XFullVector,size(X{1},1),'Replace',false);
                end 
                if i==2%length(X)
                    X=Z;break;
                end
                continue;
            end
            
            if permuteMode==8%Brute force permutation
                disp('Each condition with adding noise.')
                if i==1
                    Z{i}=X{1};%max(XFull(:,randi([1 size(XFull,2)],1,size(X{i},2))),0);  
                else
                    Z{i}=max(Z{1}+std(Z{1}')'*1*ones(1,size(Z{1},2)).*(rand(size(Z{1},1),size(Z{1},2))-0.5)*2,eps);
                    assignin('base', 'stda', std(Z{1}')');%dd
                    if i==2
                    X=Z;
                    break;
                    end
                end
                continue;
            end
            
            if permuteMode<=3
            if permuteMode==2
                if i==1
             [E,index] = sortrows(X{i},1,'descend');
                end
            elseif permuteMode==1
                if i==1
             [E,index] = sortrows(X{i},1,'ascend');
                end
            elseif permuteMode==3
                if i==1
                index=randperm(size(X{i},1));
                disp('random permute')
                %disp(index);
                end
            end
            X{i}=X{i}(index,:);
            end
            
            
        end
        %Y=X;
        
        %assignin('base', 'Y', Y);dd
        disp('finish permuation data');
       
end