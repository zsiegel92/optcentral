function a1bounds = xlsBounds( row1,col1,row2,col2 )
 a1bounds = [ExcelA1Gen(row1,col1), ':', ExcelA1Gen(row2,col2)];
end

