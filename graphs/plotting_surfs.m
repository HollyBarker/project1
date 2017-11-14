dataLARGER=load('xfinalLARGER.txt');
dataLARGERBANDED=load('xfinalLARGERBANDED.txt');

surf(dataLARGER)
zlabel('X','FontSize', 14)
title('Final x values for the CG method and MMatrix','FontSize', 14)

figure;
surf(dataLARGERBANDED)
zlabel('X','FontSize', 14)
title('Final x values for the CG method and MBandedMatrix','FontSize', 14)
