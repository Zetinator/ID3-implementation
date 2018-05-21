function [tree] = plotTree(obj,tree, auxCounter)
  auxCounter = auxCounter + 1;
  if obj.class > 0
  else
    if size(obj.descendant,2) > 0
      for i=1:size(obj.descendant,2)
        tree = [tree, auxCounter];
        % tree = [tree, plotTree(obj.descendant(i), tree, auxCounter)];
        tree = plotTree(obj.descendant(i), tree, auxCounter);

      end

    end

  end
  display('SUCCESS')
end
