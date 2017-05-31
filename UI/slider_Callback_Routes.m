function slider_Callback_Routes(hObject, eventdata, handles)
% hObject    handle to slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine...

%% Slider Value
slider_value = get(hObject,'Value');
slider_value = floor(slider_value);
set(hObject,'Value',slider_value);
%map_update(slider_value,get(hObject,'UserData'));

ud = get(hObject,'UserData');







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
%G = ud.G;
















%% Draw Graph
dx = .1; dy = .1;
p=plot(ud.G,'XData',ud.coordinates(2,:),'YData',ud.coordinates(1,:),'MarkerSize',7);

%% Edge Information: Transportation

% edgelab = cellfun(@(x) {[num2str(x),' Trucks']}, num2cell(G.Edges.Weight));
% edgeArray = cellfun(@(x) findnode(G,x),G.Edges.EndNodes);%nEdges by 2 containing edge endpoints
% edgeXMids = (arrayfun(@(x) ud.coordinates(2,x), edgeArray(:,1)) +arrayfun(@(x) ud.coordinates(2,x), edgeArray(:,2)))/2;
% edgeYMids = (arrayfun(@(x) ud.coordinates(1,x), edgeArray(:,1)) +arrayfun(@(x) ud.coordinates(1,x), edgeArray(:,2)))/2;

% if (size(ud.edgelab{slider_value},1)>0)
%     %text(ud.edgeXMids{slider_value}+dx,ud.edgeYMids{slider_value}+dy,ud.edgelab{slider_value});
% end
text(ud.edgeXMids+dx,ud.edgeYMids+dy,ud.edgelab);

%% Labeling node information: production, sales, (holding?)
% Nothing yet...

%% Coloring Nodes
%colormap parula;

%% Editing text box

for i = 1:length(hObject.Parent.Children)
    if isprop(hObject.Parent.Children(i),'Style')
        if strcmp(hObject.Parent.Children(i).Style,'text')
            if slider_value>0
                hObject.Parent.Children(i).String = ['Time-step: ',num2str(slider_value)];
            else
                hObject.Parent.Children(i).String = ['All-time Statistics'];
            end
        end
    end
end

end
