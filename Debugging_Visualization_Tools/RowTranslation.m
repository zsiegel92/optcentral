function [] = RowTranslation(constraintRow,constraintMatrix,Jinv)
variableInds = find(constraintMatrix(constraintRow,:));
variables = Jinv(variableInds,:);
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
variableNames = {'y','z','t','w','z0','t0','z00','wm','wm0'};
variableEntry = [];
for i = 1:size(variables,1)
    variableEntry = strcat(num2str(constraintMatrix(constraintRow,variableInds(i))),'*',variableNames(variables(i,1)),'(',num2str(variables(i,2)),',',num2str(variables(i,3)),',',num2str(variables(i,4)),',',num2str(variables(i,5)),')');
    disp(char(variableEntry))
end

end