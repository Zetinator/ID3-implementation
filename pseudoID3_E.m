function [descendant] = pseudoID3_E(data,featureLabels,ID,featureToCheck,classMark,width)
  display(data)
  % ID3 code
  % The class always in the first column... the rest of the columns are features...
  features = data(:,2:size(data,2));
  labels = data(:,1);

  % terminal case
  if isempty(data)
  else
    numOfClass = unique(labels);
    if size(numOfClass,1) == 1
      % leaf founded
      strID = sprintf('%s', dicomuid);
      descendant = [id3(ID, strID, featureToCheck, classMark, width, numOfClass)];
      fprintf('node: %s found class: %d \n', descendant(1).ID, numOfClass);
    else
      % Sturges rule to select a good number of intervals for the histogram.
      numInterval = ceil(1+log2(size(nonzeros(features(:,1)),1)));
      amplitud = range(features)/(numInterval);

      % Number of features of value j in the range of the interval j;
      for i=1:size(features,2)
        [counts,centers] = hist(features(:,i),numInterval);
        % Should be an interval between  and numInterval??? ##Debug later
        occurrences(i,:) = counts;
        classMarks(:,i) = centers';
      end

      % Transpose and normalize to compute probability
      probability = (1/size(features,1)) * occurrences';

      % Data per class (in submatrix)
      for i = 1:size(data,1)
        for j = 1:size(numOfClass,1)
          if data(i,1) == numOfClass(j)
            % Maybe subMat(j,i,:) its more intuitive
            subMat(i,:,j) = data(i,2:size(data,2));
          end
        end
      end

      % Probability of class Cj when feature Ak has a value of j
      for j = 1:size(numOfClass,1)
        for i=1:size(subMat,2)
          occurrencesPerClass(:,i,j) = hist(nonzeros(subMat(:,i,j)),classMarks(:,i))';
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

      % Label of father node would be the position index of the feature with min entropy...
      featureToCheck = featureLabels(pos);

      % prepare the new featureLabels that are going to be passed to the recursive children of the function
      cutLL = featureLabels(1:pos-1);
      cutLR = featureLabels((pos+1):end);
      featureLabels = [cutLL,cutLR];

      % New data that is going to be passed to the recursive children of the function...
      cutIzq = data(:,1:pos);
      cutDer = data(:,(pos+2):end);
      data = [cutIzq,cutDer];

      % For each attribute of the smaller entropy I take submatrices of dataFeatures for the "j" intervals to feed the "j" children
      for i = 1:size(data,1)
        for j = 1:numInterval
          % relative epsilon to avoid floating point aritmetics mismatches
          epsilon = abs(features(i,pos) - classMarks(j,pos))/100;
          % check if the data belongs to the interval "j" using the distance to the class mark and range as threshold
          if abs(features(i,pos) - classMarks(j,pos)) < (amplitud(pos)/2 + epsilon)
            subMatrixs(i,:,j) = data(i,:);
          end
        end
      end

      strID = sprintf('%s', dicomuid);
      node = id3(ID, strID, featureToCheck, classMark, width, 0);
      for i = 1:numInterval
        clear bringIt
        auxCounter = 1;

        % Clean the matrix of the zeros... can be debug later...
        for k = 1:size(features,1)
          if any(subMatrixs(k,:,i)')
            bringIt(auxCounter,:) = subMatrixs(k,:,i);
            auxCounter = auxCounter + 1;
          end

        end
        if any(subMatrixs(:,:,i))
          classMark = classMarks(i,pos);
          width = amplitud(pos)/2;
          % pushback in the descendant vector
          node.descendant = [node.descendant, pseudoID3_E(bringIt,featureLabels,strID, featureToCheck,classMark,width)];
        end
      end
      descendant = node;

    end
  end
end
