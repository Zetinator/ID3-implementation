% MasterID3
function [ trees ] = MasterID3(data)
	labels = [1:size(data,2)-1];
	trees = pseudoID3(data,labels,0,0,0,0,0,[]);
	display(vpa(trees,5));

	% treeplot
	for i = 1: size(trees,1)
		nodos(i) = trees(i,1) - min(trees(:,1));
	end
	nodos = sort(nodos);
	display(nodos');
	treeplot(nodos);
end

