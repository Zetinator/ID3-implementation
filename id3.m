classdef id3 < handle
  properties
    fatherID;
    ID;
    featureToCheck;
    classMark;
    width;

    descendant = [];

    class;


  end


  methods
    % constructor
    function obj = id3(fatherID, ID, featureToCheck, classMark, width, class)
        if nargin == 0
          obj.fatherID = 0;
          obj.ID = 0;
          obj.featureToCheck = 0;
          obj.classMark = 0;
          obj.width = 0;
          obj.class = 0;

        elseif nargin == 6
          obj.fatherID = fatherID;
          obj.ID = ID;
          obj.featureToCheck = featureToCheck;
          obj.classMark = classMark;
          obj.width = width;
          obj.class = class;

        else
          error('Wrong number of input arguments')
        end
    end

    function [success] = createTree(obj, data)
      % When given the data matrix the "class" always in the first column... the rest of the columns are features...
      featureLabels = [1:size(data,2)-1];
      obj.descendant = [];
      obj.descendant = [obj.descendant, pseudoID3_E(data,featureLabels,0,0,0,0)]

      display('SUCCESS')
    end

    function [success] = plot(obj)
      % plot the tree model graph
      tree = plotTree(obj,[],0);
      tree = tree - 1;
      treeplot(tree);

      display('SUCCESS')
    end

  end


end
