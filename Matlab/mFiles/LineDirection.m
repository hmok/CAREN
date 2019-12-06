a = [1 3 4 1 -1.5 2 -3 .5];% e.g. choose your Y direction
b= [-1 2 -3 -1 1 -3 1 .5];% choose your X direction
for i = 1:length(a)

    a1 =[i i+a(i)];
    b2 = [0 b(i)];

    line(a1,b2) 
    hold on
end