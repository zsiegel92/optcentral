function popup_Callback(hObject, eventdata, handles)
ud = get(hObject,'UserData');
popupVal = hObject.Value;
%options = hObject.String;


if (popupVal == 1)
    cla;
    axesHandlesToChildObjects = findobj('Tag', '50States');
    delete(axesHandlesToChildObjects);
    
    
    h = findobj('Tag', 'slider_zt');
    delete(h);
    h = findobj('Tag', 'slider_Routes');
    delete(h);
    h = findobj('Tag','text_zt');
    delete(h);
    h=findobj('Tag','text_Splash');
    delete(h);
    uicontrol('Style','text',...
        'Position',[20 100 500 200],...
        'Tag','text_Splash',...
        'String','Choose a functionality from the drop-down menu.');
    %     h = findobj('Tag','50States');
    %     delete(h)
    
    
elseif (popupVal ==2)
    
    h = findobj('Tag', 'slider_zt');
    delete(h);
    h = findobj('Tag','text_zt');
    delete(h);
    h=findobj('Tag','text_Splash');
    delete(h);
    h = findobj('Tag', 'slider_Routes');
    delete(h);
    
    Transport_and_Production_Graphing_Animation;
    %     f = figure('Visible','off');
    %     popUpMenu = uicontrol('Style', 'popupmenu',...
    %         'String', {'Transport and Production','Route Planning and Inventory','Production Capacity','Route and Production Statistics Map'},...
    %         'Position', [20 340 100 50],...
    %         'Callback', @popup_Callback,...
    %         'UserData',ud);
    %     set(popUpMenu,'Value',popupVal);
    %
    %     f.Visible = 'on';
    
elseif (popupVal ==3)
    h = findobj('Tag', 'slider_zt');
    delete(h);
    h = findobj('Tag','text_zt');
    delete(h);
    h=findobj('Tag','text_Splash');
    delete(h);
    h = findobj('Tag', 'slider_Routes');
    delete(h);
    Route_Planning;
    
    %     f.Visible = 'on';
elseif (popupVal ==4)
    h = findobj('Tag', 'slider_zt');
    delete(h);
    h = findobj('Tag','text_zt');
    delete(h);
    h=findobj('Tag','text_Splash');
    delete(h);
    h = findobj('Tag', 'slider_Routes');
    delete(h);
    clear hi;
    
    BOMVisualizer(ud.Gbom);
end
end

