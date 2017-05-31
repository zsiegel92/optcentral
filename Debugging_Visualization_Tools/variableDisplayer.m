%{
%Index Description
Populate Index Matrix J such that:
J(1,i,j,k,1) = io y(i,j,k) = quantity of i held at location k at time j
J(2,i,j,k,1) = io z(i,j,k) = quantity of i produced at location k at time j
J(3,i,j,k,L) = io t(i,j,k,L) = quantity of i transported from location k  to L at time j
J(4,i,j,k,1) = io w(i,j,k) = quantity of i sold at location k at time j
J(5,i,j,k,1) = io z0(i,j,k) = (1 or 0) whether i is produced at location k
at time j
J(6,1,j,k,L) = t0(j,k,L) = number of trucks from location k to L at time
j
J(7,i,j,k,1) = io z00(i,j,k) = (1 or 0) whether i is "in production"
(at least one time-step after start)at location k and time j
J(8,i,j,k,m)= wm(i,j,k)
J(9,i,j,k,m)=wm0(i,j,k)
(io means "index of")
%}

[JT,JinvT,intconT]= Indexer(param);
[qIneqT,bIneqT,qEqT,bEqT,qObT,LBT,UBT] = constraint_formalizer(param,JT);
% 
% JinvT = Jinv2;
qT = qEqT;
% JT = J2;
VariableInQuestion = [2,5,3,1,1];

variableNames = {'y','z','t','w','z0','t0','z00','wm','wm0'};
VIQ = VariableInQuestion;
disp(char(strcat({'Rows Containing '}, variableNames(VariableInQuestion(1)),'(',num2str(VariableInQuestion(2)),',',num2str(VariableInQuestion(3)),',',num2str(VariableInQuestion(4)),',',num2str(VariableInQuestion(5)),'):')));
indexOfVariable = JT(VIQ(1),VIQ(2),VIQ(3),VIQ(4),VIQ(5));
inds= find(qT(:,indexOfVariable)); %returns [355, 357]
for i = 1:length(inds)
    disp(['Row number ',num2str(inds(i)),':']);
    %JinvT(find(qT(inds(1),:)),:)
    RowTranslation(inds(i),qT,JinvT)
    disp(sprintf ('\n') );
end

% prunerT = Pruning(param,JT);
% nnz(prunerT==indexOfVariable)


% prunedVars =JinvT(prunerT,:);