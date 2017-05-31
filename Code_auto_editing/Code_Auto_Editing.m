%Trying to replace parameters of the form 'PARAMETER' with parameters of the
%form 'param.PARAMETER'


%param = Manual_Parameterizer();

names = fieldnames(param);
newNames = cellfun(@(x) ['param.',x],names,'UniformOutput',0);
cd ~/Desktop

fileID = fopen('OptCentral_tiered_TEXT.txt'); %,'%s','Delimiter','\n'
myCode = textscan(fileID,'%s','Delimiter','\n');
myCode=myCode{1,1};

for i = 4:length(names) %omit [nz; nt; nloc] = names(1:3)
    if (~any(i == [1,2,3,35,31]))
        myCode = strrep(myCode,char(names(i)),['param.',char(names(i))]);
    end
end
fclose(fileID);
fileID = fopen('Code_auto_editing/OptCentral_tiered_TEXT_EDITED.txt','w');
formatSpec = '%s\n';
[nrows,ncols]=size(myCode);
for row = 1:nrows
    fprintf(fileID,formatSpec,myCode{row,:});
end
fclose(fileID);