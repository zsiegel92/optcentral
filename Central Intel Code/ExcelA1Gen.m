function a1name = ExcelA1Gen( rows,cols )
 a1name = [char(base2dec(dec2base(cols,26)',26)+64)' num2str(rows)];
end

