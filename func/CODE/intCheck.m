function intCheck(in)

integer = sum(abs(in-round(in)) > 1e-8) == 0;
 
if integer == 0  
    error('inputs are not integer valued, returning')
else
    display('inputs are integer valued, you can proceed')
end