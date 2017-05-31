function slider_Callback_ZT(hObject, eventdata, handles)
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



%% Draw Graph
dx = .1; dy = .1;
% G = digraph(squeeze(ud.solution.T0(slider_value,:,:)),ud.locNames); %Add number of trucks travelling each route
%p=plot(G,'XData',coordinates(2,:),'YData',coordinates(1,:),'EdgeCData',zeros(1,nEdges));
% 
subplot(2,1,1);
p = plot(ud.G{slider_value},'XData',ud.coordinates(2,:),'YData',ud.coordinates(1,:),'MarkerSize',7);


%'Time-Step','Origin','Destination','Number Trucks Sent'	
subplot(2,1,2)
% t = uitable(hObject.Parent,'Data',ud.solution.track_T0,'ColumnName',{'Time-Step','Origin','Destination','Number Trucks Sent'},'Position',[20 20 400 100]);

%% Edge Information: Transportation

% edgelab = cellfun(@(x) {[num2str(x),' Trucks']}, num2cell(G.Edges.Weight));
% edgeArray = cellfun(@(x) findnode(G,x),G.Edges.EndNodes);%nEdges by 2 containing edge endpoints
% edgeXMids = (arrayfun(@(x) ud.coordinates(2,x), edgeArray(:,1)) +arrayfun(@(x) ud.coordinates(2,x), edgeArray(:,2)))/2;
% edgeYMids = (arrayfun(@(x) ud.coordinates(1,x), edgeArray(:,1)) +arrayfun(@(x) ud.coordinates(1,x), edgeArray(:,2)))/2;

if (size(ud.edgelab{slider_value},1)>0)
    text(ud.edgeXMids{slider_value}+dx,ud.edgeYMids{slider_value}+dy,ud.edgelab{slider_value});
end


%% Labeling node information: production, sales, (holding?)

dx = .25; dy = .4;
text(ud.coordinates(2,ud.zLocs{slider_value})+dx,ud.coordinates(1,ud.zLocs{slider_value}) + dy, ud.nodelab{slider_value});


%% Coloring Nodes
%colormap parula;

caxis([0 1]);
colorbar('FontWeight','bold','FontSize',12,'Ticks',[0,0.5,1],'TickLabels',{'Idle','Initiating Production','In-Production'});

% nodecolors =ud.solution.Z00(slider_value,:) + 0.5*sum(squeeze(ud.solution.Z0(:,slider_value,:)),1); %Is anything starting production?
p.NodeCData = ud.nodecolors{slider_value};

%% Editing text box

for i = 1:length(hObject.Parent.Children)
    if isprop(hObject.Parent.Children(i),'Style')
        if strcmp(hObject.Parent.Children(i).Style,'text')
            hObject.Parent.Children(i).String = ['Time-step: ',num2str(slider_value)];
        end
    end
end

end
