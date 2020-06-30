function X=PermDataInputOrder(X,permuteMode,index)

        assignin('base', 'X1', X);
        for i=1:length(X)    
            X{i}=X{i}(index,:);
        end
        disp('finish permuation data by inputting orders');
       
end