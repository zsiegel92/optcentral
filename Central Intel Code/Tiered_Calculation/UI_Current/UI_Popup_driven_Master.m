% Before running this script, either run OptCentral_tiered, or load a
% workspace saved after running OptCentral.
function [] = UI_Popup_driven_Master(oldSol)
nt = oldSol.params.nt;
nloc = oldSol.params.nloc;
nz = oldSol.params.nz;
%% Initializing figure (w/ US outline in background)
f = figure('Visible','off');

%% Creating popupControl uiControl
userData_popup =struct('nt',nt,'nloc',nloc,'coordinates',oldSol.param.coordinates);%All variables necessary for updating plot (T, Z, W, etc.)
% %Have to add cell array parameters separately because cell arrays behave differently in structure constructor.
userData_popup.locNames=oldSol.params.locNames;
userData_popup.solution = oldSol.SOLUTION;
userData_popup.params = oldSol.params;
userData_popup.Gbom = oldSol.params.Gbom;
popUpMenu = uicontrol('Style', 'popupmenu',...
    'String', {'Blank','Transport and Production','Route Planning and Inventory','Bill of Materials (BOM)'},... %,'Production Capacity','Route and Production Statistics Map'},...
    'Position', [20 340 100 50],...
    'Callback', @popup_Callback,...
    'UserData',userData_popup);




%% Initializing figure content for dropdown menu state 1
% ud = userData_popup;
% Transport_and_Production_Graphing_Animation;
txt=uicontrol('Style','text',...
    'Position',[20 100 500 200],...
    'Tag','text_Splash',...
    'String','Choose a functionality from the drop-down menu.');




%% Making graphic visible

f.Visible = 'on';

end