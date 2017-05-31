function [BOM,Gbom] = BOMCREATOR()
%Uses exactly 11 products
%Call with:
%[p,BOM,Gbom]=BOMCREATOR;
Gbom = digraph;
Gbom = addnode(Gbom,{'Wheat Berries (lb)', 'Beef (lb)', 'Egg (ea)','Salt (g)','Wheat Flour (lb)','Ground Beef (lb)','Roll/Bun (ea)', 'Burger (8) FG','Egg Sandwich (ea) FG'});
%disp([array2table((1:numnodes(Gbom))'), Gbom.Nodes]);
Gbom = addedge(Gbom,{'Wheat Berries (lb)'},{'Wheat Flour (lb)'},[1]);

Gbom = addedge(Gbom,table([5,7;4,7;4,8;6,8;2,6;7,8;7,9;3,9],[0.1,1,0.1,.25,1,1,1,1]','VariableNames',{'EndNodes','Weight'}));


%Adding More Nodes and edges
Gbom = addnode(Gbom,'Mayonnaise (lb)');
Gbom = addnode(Gbom,'Vegetable Oil (lb)');
Gbom = addedge(Gbom,table([3,10;11,10;10,9;10,8;11,9],[2,3/4,1/16,1/16,1/32]','VariableNames',{'EndNodes','Weight'}));


%Exporting to BOM matrix
[s,t]=findedge(Gbom);
BOM = full(sparse(s,t,Gbom.Edges.Weight,numnodes(Gbom),numnodes(Gbom)))';
clear s; clear t;


%Defining Raw Materials and Finished Goods
Gbom.Nodes.rawMaterials = zeros(11,1);
Gbom.Nodes.FG=zeros(11,1);%SKUS that are raw materials

%Manually label raw materials and FG
%Gbom.Nodes.rawMaterials([1,2,3,4,11])=1
%Gbom.Nodes.FG([8,9])=1;%Skus that are finished goods
for i = 1:size(BOM,1)
    if (nnz(BOM(i,:))==0)
        Gbom.Nodes.rawMaterials(i)=1; %Automatically label raw materials
    end
    if (nnz(BOM(:,i))==0)
        Gbom.Nodes.FG(i)=1; %Automatically label finished goods
    end
end

%Manually label finished goods that can also be sub-assemblies
Gbom.Nodes.FG(7)=1; %Buns/Rolls

%disp([array2table((1:numnodes(Gbom))'), Gbom.Nodes]); %List node names and
                                                        %properties

%Visualize Gbom
%p=BOMVisualizer(Gbom);

end