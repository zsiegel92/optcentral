function p = BOMVisualizer(Gbom)

%Possibly a way to clear or replace background of image plot.
% ha = axes('units','normalized','position',[0,0,1,1]);
% % Move the background axes to the bottom
% uistack(ha,'bottom');
% % Load in a background image and display it using the correct colors
% hi = imagesc(imread('BO.jpg'),'Tag','Bo');
% colormap gray
% % Turn the handlevisibility off so that we don't inadvertently plot into the axes again
% % Also, make the axes invisible
% ha.GridAlpha=0;
% ha.MinorGridAlpha = 0;
% set(ha,'handlevisibility','off', ...
%     'visible','off');
% 


p=plot(Gbom,'EdgeLabel',Gbom.Edges.Weight,'NodeLabel',Gbom.Nodes.Name);


highlight(p,find(Gbom.Nodes.rawMaterials),'NodeColor','g'); %Highlight raw ingredients in GREEN
highlight(p,find(Gbom.Nodes.FG),'NodeColor','r'); %Highlight FG in RED


% AUTOMATE THIS USING p.XData and p.Ydata
LegendX = 5.5;
LegendY = 4.25;
dx = .2;
dy = .25;

nodeLegend = {{'Raw Material','\bullet','green'},{'Finished Good','\bullet','red'},{'Sub-Assembly','\bullet','blue'},{sprintf('(Amount of) \nIngredient Used'),'\rightarrow','blue'}};

for i = 1:length(nodeLegend)
    text(LegendX, LegendY -(i-1)*dy, nodeLegend{i}{2}, 'Color', nodeLegend{i}{3});
    text(LegendX +dx, LegendY - (i-1)*dy, nodeLegend{i}{1}, 'Color', 'black');
end

title('Bill Of Materials')
end