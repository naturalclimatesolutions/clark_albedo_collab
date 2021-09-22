function LeapYear = isLeapYear(intYear)
a = mod(intYear,400);
b = mod(intYear,100);
c = mod(intYear,4);
LeapYear = ((a==0)||(b~=0 && c==0));
end