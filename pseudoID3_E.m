function [tree] = pseudoID3(data,labels,nivel,labelFather,branch,brandOfTheCurse,greatness,tree)
% ID3
% In which level am I, of the recursive tree?
nivel = nivel+1;
% The class always in the first column...
features = data(:,2:size(data,2));
theClass = data(:,1);

% Probability of attribute Ak to have a value J
if isempty(data)
else
    numOfClass = unique(theClass);
    if size(numOfClass,1) == 1
        tree(end+1,:) = [nivel,labelFather,branch,brandOfTheCurse,greatness,numOfClass]; 
        
    else
        % Sturges rule
        numInterval = ceil(1+log2(size(nonzeros(features(:,1)),1)));
        amplitud = range(features)/(numInterval);
        
        % Number of features of value j in the range of the interval j;
        for i=1:size(features,2)
            [counts,centers] = hist(features(:,i),numInterval);
            occurrences(i,:) = counts;
            classMarks(:,i) = centers';
        end
        
        % Transpose and normalize
        probability = (1/size(features,1)) * occurrences';
        
        % Data per class (in submatrix)
        for i = 1:size(data,1)
            for j = 1:size(numOfClass,1)
                if data(i,1) == numOfClass(j)
                    subMat(i,:,j) = data(i,:);
                end
            end
        end
        
        subM = subMat(:,2:size(subMat,2),:);

        % Probability of class Cj when Ak has a value of j
        for j = 1:size(numOfClass,1)
            for i=1:size(subM,2)
                occurrencesPerClass(:,i,j) = hist(nonzeros(subM(:,i,j)),classMarks(:,i))';
            end
            probabilityPerClass(:,:,j) = occurrencesPerClass(:,:,j)/sum(occurrencesPerClass(:,1,j));
        end
        
        % Entrophy of feature Ak, by Claude Shannon
        for j = 1:size(features,2)
            dEntrophy = 0;
            for i = 1:numInterval
                alfa = 0;
                for k = 1:size(numOfClass,1)
                    if probabilityPerClass(i,j,k) == 0
                        aux = 1;
                    else
                        aux = probabilityPerClass(i,j,k);
                    end
                    alfa = probabilityPerClass(i,j,k) * log2(aux) + alfa;
                end
                dEntrophy = dEntrophy + probability(i,j) * (-1) * alfa;
            end
            entrophy(j) = dEntrophy;
        end
        

        % Arrangement against case where the entropy is 0 is taken the attribute of the right...
        if nnz(entrophy) == 0
            pos = 1;

        else
            [minimum,pos] = min(entrophy);
        end

        % Label of father node...
        labelFather = labels(pos);
        cutLL = labels(1:pos-1);
        cutLR = labels((pos+1):end);
        labels = [cutLL,cutLR];

        % New data...
        cutIzq = data(:,1:pos);
        cutDer = data(:,(pos+2):end);
        data = [cutIzq,cutDer];
        
        % For each attribute of smaller entropy I take submatrices for the j ranges
        for i = 1:size(data,1)
            for j = 1:numInterval
                epsilon = abs(features(i,pos) - classMarks(j,pos))/100;
                if abs(features(i,pos) - classMarks(j,pos)) < (amplitud(pos)/2 + epsilon)
                    subMatrixs(i,:,j) = data(i,:);
                end
            end
        end
        
        % RecursivityTime
        for i = 1:numInterval
            clear bringIt
            for k = 1:size(features,2)
                bringIt(:,k) = nonzeros(subMatrixs(:,k,i));
            end
            branch = i;
            brandOfTheCurse = classMarks(i,pos);
            greatness = amplitud(pos);
            tree(end+1,:) = [nivel,labelFather,branch,brandOfTheCurse,greatness,0]; 
            tree = pseudoID3(bringIt,labels,nivel,labelFather,branch,brandOfTheCurse,greatness,tree);
        end
        
    end
end
end

