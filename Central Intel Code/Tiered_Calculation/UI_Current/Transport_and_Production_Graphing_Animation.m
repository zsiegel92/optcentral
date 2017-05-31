%% This script is called from 'popup_Callback_ZT', so only have access to information in ud = get(hObject,'UserData')


%% Getting necessary userData elements.
nt = ud.nt;
coordinates = ud.coordinates;

%% Creating the 'background' figure and axes
title('Transportation and Production Plan')
ud.ha = axes('units','normalized','position',[0,0,1,1]);
% Move the background axes to the bottom
uistack(ud.ha,'bottom');
% Load in a background image and display it using the correct colors
% The image used below, is in the Image Processing Toolbox.  If you do not have %access to this toolbox, you can use another image file instead.
ud.hi = imagesc(imread('US_outline.jpeg'),'Tag','50States');
colormap gray
% uistack(hi,'bottom');
% Turn the handlevisibility off so that we don't inadvertently plot into the axes again
% Also, make the axes invisible
ud.ha.GridAlpha=0;
ud.ha.MinorGridAlpha = 0;
set(ud.ha,'handlevisibility','off', ...
            'visible','off');

%% Initializing graph; Generating edge labels; Finding midpoints of non-zero edges
%G = digraph(squeeze(ud.solution.T0(1,:,:)),locNames);
G = cell(nt,1);
edgelab = cell(nt,1);
edgeXMids = cell(nt,1);
edgeYMids = cell(nt,1);
for i = 1:nt
    G{i} = digraph(squeeze(ud.solution.T0(i,:,:)),ud.params.locNames);
    edgelab{i} = cellfun(@(x) {[num2str(x),' Trucks']}, num2cell(G{i}.Edges.Weight));
    
    edgeArray = cellfun(@(x) findnode(G{i},x),G{i}.Edges.EndNodes); %nEdges by 2 containing edge endpoints
    edgeXMids{i} = (arrayfun(@(x) coordinates(2,x), edgeArray(:,1)) +arrayfun(@(x) coordinates(2,x), edgeArray(:,2)))/2;
    edgeYMids{i} = (arrayfun(@(x) coordinates(1,x), edgeArray(:,1)) +arrayfun(@(x) coordinates(1,x), edgeArray(:,2)))/2;
end

%% Adding text to edges
dx = .1; dy=.1;

if (size(edgelab{1},1)>0)
    text(edgeXMids{1}+dx,edgeYMids{1}+dy,edgelab{1});
end

%% Labeling node information: production, sales, (holding?)

dx = .25; dy = .4;
zProducts = cell(nt,1);
zLocs = cell(nt,1);
nodelab = cell(nt,1);

for i = 1:nt
[zProducts{i},zLocs{i}] = find(squeeze(ud.solution.Z0(:,i,:)));
% + ud.solution.Z00(:,slider_value,:)
nodelab{i} = cellfun(@(x) ['Producing sku: ',num2str(x)],num2cell(zProducts{i}),'UniformOutput',0);



end


text(coordinates(2,zLocs{1})+dx,coordinates(1,zLocs{1}) + dy, nodelab{i});






%% Plotting graph




subplot(2,1,1);
p=plot(G{1},'XData',coordinates(2,:),'YData',coordinates(1,:),'MarkerSize',7);

% info = ud.solution.track_T0;
info = ud.solution.track_T;
yPerRow = 15;
tableHeight = size(info,1)*yPerRow+yPerRow;
tableWidth = 374; %for 4 columns...
ColumnName = {'Time-Step','Origin','Destination','SKU','Amount Transported'};
xPerChar = 10;
tableWidth=length(cell2mat(ColumnName))*xPerChar;


subplot(2,1,2)
t = uitable(hObject.Parent,'Data',info);
%t.ColumnName={'Time-Step','Origin','Destination','Number Trucks Sent'};
t.ColumnName=ColumnName;
t.Position= [40; 30+min(max(0,tableHeight-40),40); tableWidth; tableHeight];
axis off
%% Coloring Nodes

nodecolors = cell(nt,1);
for i = 1:nt
nodecolors{i} = ud.solution.Z00(i,:) + 0.5*sum(squeeze(ud.solution.Z0(:,i,:)),1); %Is anything starting production?
end

colormap parula;

p.NodeCData = nodecolors{1};

caxis([0 1]);
colorbar('Ticks',[0,0.5,1],'TickLabels',{'Idle','Initiating Production','In-Production'},'FontWeight','bold','FontSize',12);



%% Creating Slider uiControl
userData_slider_ZT =struct('coordinates',coordinates);%All variables necessary for updating plot (T, Z, W, etc.)
%Have to add cell array parameters separately because cell arrays behave differently in structure constructor.
userData_slider_ZT.locNames=ud.params.locNames;
userData_slider_ZT.solution = ud.solution;
userData_slider_ZT.G = G;
userData_slider_ZT.edgelab = edgelab;
userData_slider_ZT.edgeXMids = edgeXMids;
userData_slider_ZT.edgeYMids = edgeYMids;
userData_slider_ZT.nodelab = nodelab;
userData_slider_ZT.zLocs = zLocs;
userData_slider_ZT.nodecolors = nodecolors;
sld = uicontrol('Style', 'slider',...
        'Min',1,'Max',nt,'Value',1,...
        'Position', [400 20 120 20],...
        'Callback', @slider_Callback_ZT,...
        'Tag','slider_zt',...
        'UserData',userData_slider_ZT); 

set(sld,'SliderStep',[1/(ud.params.nt-1),1/(ud.params.nt-1)]) %To improve
%slider behavior

%% Creating time-step 'edit' textbox uiControl
txt=uicontrol('Style','text',...
       'Position',[400 45 120 20],...
       'Tag','text_zt',...
        'String',['Time-step:',num2str(1)]);

