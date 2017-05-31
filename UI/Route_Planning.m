%% This script is called from 'popup_Callback_ZT', so only have access to information in ud = get(hObject,'UserData') (from popup_Callback)


%% Getting necessary userData elements.

nt = ud.nt;
coordinates = ud.coordinates;

% This creates the 'background' axes
ha = axes('units','normalized','position',[0,0,1,1]);
% Move the background axes to the bottom
uistack(ha,'bottom');
% Load in a background image and display it using the correct colors
% The image used below, is in the Image Processing Toolbox.  If you do not have %access to this toolbox, you can use another image file instead.
hi = imagesc(imread('US_outline.jpeg'),'Tag','50States');
colormap gray
% Turn the handlevisibility off so that we don't inadvertently plot into the axes again
% Also, make the axes invisible
ha.GridAlpha=0;
ha.MinorGridAlpha = 0;
set(ha,'handlevisibility','off', ...
    'visible','off');
title('Route Planning and Inventory')
%% Initializing graph; Generating edge labels; Finding midpoints of non-zero edges
%G = digraph(squeeze(ud.solution.T0(1,:,:)),locNames);

[s,t]=find(squeeze(sum(ud.solution.T0,1)));%The routes that are actually travelled within the system.
edgeInfo = cell(length(s),1); %i'th row consists of [R_t(s(i),t(i)),f_t(s(i),t(i))]
for i=1:length(s)
    edgeInfo{i}=[ud.params.R_t(s(i),t(i)),ud.params.f_t(s(i),t(i))];
end

G = digraph(s,t,zeros(length(s),1),ud.params.locNames);
%edgelab = cellfun(@(x) {[num2str(x),' Trucks']}, num2cell(G.Edges.Weight));

% edgeArray = cellfun(@(x) findnode(G,x),G.Edges.EndNodes); %nEdges by 2 containing edge endpoints
% edgeXMids = (arrayfun(@(x) coordinates(2,x), edgeArray(:,1)) +arrayfun(@(x) coordinates(2,x), edgeArray(:,2)))/2;
% edgeYMids = (arrayfun(@(x) coordinates(1,x), edgeArray(:,1)) +arrayfun(@(x) coordinates(1,x), edgeArray(:,2)))/2;


% edgeArray = cellfun(@(x) findnode(G,x),G.Edges.EndNodes); %nEdges by 2 containing edge endpoints
edgeXMids = (arrayfun(@(x) coordinates(2,x), s) +arrayfun(@(x) coordinates(2,x),t))/2;
edgeYMids = (arrayfun(@(x) coordinates(1,x), s) +arrayfun(@(x) coordinates(1,x),t))/2;


%edgelab = cellfun(@(x)['Route time: ', 
edgelab = cellfun(@(x) {sprintf(['R_t: ',num2str(x(1)),'\nf_t: ',num2str(x(2))])},edgeInfo);



%% Labeling node information: production, sales, (holding?)
% 
% dx = .25; dy = .4;
% zProducts = cell(nt,1);
% zLocs = cell(nt,1);
% nodelab = cell(nt,1);
% 
% for i = 1:nt
% [zProducts{i},zLocs{i}] = find(squeeze(ud.solution.Z0(:,i,:)));
% % + ud.solution.Z00(:,slider_value,:)
% nodelab{i} = cellfun(@(x) ['Producing sku: ',num2str(x)],num2cell(zProducts{i}),'UniformOutput',0);
% 
% 
% 
% end
% 
% 
% text(coordinates(2,zLocs{1})+dx,coordinates(1,zLocs{1}) + dy, nodelab{i});
% 
% 




%% Plotting graph
p=plot(G,'XData',coordinates(2,:),'YData',coordinates(1,:),'MarkerSize',7);

%% Adding text to edges

dx = 4;
dy=.1;
text(edgeXMids+dx,edgeYMids+dy,edgelab);




%% Coloring Nodes
% 
% nodecolors = cell(nt,1);
% for i = 1:nt
% nodecolors{i} = ud.solution.Z00(i,:) + 0.5*sum(squeeze(ud.solution.Z0(:,i,:)),1); %Is anything starting production?
% end
% 
% colormap parula;
% 
% p.NodeCData = nodecolors{1};
% 
% caxis([0 1]);
% colorbar('Ticks',[0,0.5,1],'TickLabels',{'Idle','Initiating Production','In-Production'},'FontWeight','bold','FontSize',12);
% 


%% Creating Slider uiControl
% userData_slider_ZT =struct('coordinates',coordinates);%All variables necessary for updating plot (T, Z, W, etc.)
% %Have to add cell array parameters separately because cell arrays behave differently in structure constructor.
% userData_slider_ZT.locNames=ud.params.locNames;
% userData_slider_ZT.solution = ud.solution;
% userData_slider_ZT.G = G;
% userData_slider_ZT.edgelab = edgelab;
% userData_slider_ZT.edgeXMids = edgeXMids;
% userData_slider_ZT.edgeYMids = edgeYMids;
% userData_slider_ZT.nodelab = nodelab;
% userData_slider_ZT.zLocs = zLocs;
% userData_slider_ZT.nodecolors = nodecolors;
% sld = uicontrol('Style', 'slider',...
%         'Min',1,'Max',nt,'Value',1,...
%         'Position', [400 20 120 20],...
%         'Callback', @slider_Callback_ZT,...
%         'UserData',userData_slider_ZT); 
%    
% set(sld,'SliderStep',[1/(ud.params.nt-1),1/(ud.params.nt-1)]) %To improve
% %slider behavior

%% Creating time-step 'edit' textbox uiControl
% txt=uicontrol('Style','text',...
%        'Position',[400 45 120 20],...
%         'String',['Time-step:',num2str(1)]);
% 
% 

    
%% Creating Slider uiControl
userData_slider_Routes =struct('coordinates',coordinates);%All variables necessary for updating plot (T, Z, W, etc.)
%Have to add cell array parameters separately because cell arrays behave differently in structure constructor.
userData_slider_Routes.locNames=ud.params.locNames;
userData_slider_Routes.solution = ud.solution;
userData_slider_Routes.G = G;
userData_slider_Routes.edgelab = edgelab;
userData_slider_Routes.edgeXMids = edgeXMids;
userData_slider_Routes.edgeYMids = edgeYMids;
userData_slider_Routes.nt = ud.nt;
%userData_slider_Routes.nodelab = nodelab;
%userData_slider_Routes.zLocs = zLocs;
%userData_slider_Routes.nodecolors = nodecolors;
sld = uicontrol('Style', 'slider',...
        'Min',0,'Max',nt,'Value',1,...
        'Position', [400 20 120 20],...
        'Callback', @slider_Callback_Routes,...
        'Tag','slider_Routes',...
        'UserData',userData_slider_Routes); 

set(sld,'SliderStep',[1/(ud.params.nt),1/(ud.params.nt)]) %To improve
%slider behavior

%% Creating time-step 'edit' textbox uiControl
txt=uicontrol('Style','text',...
       'Position',[400 45 120 20],...
       'Tag','text_Routes',...
        'String',['Time-step:',num2str(1)]);

