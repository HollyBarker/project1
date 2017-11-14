vector10=[0.0413223,0.0743802,0.0991736,0.115702,0.123967,0.123967,0.115702,0.0991736,0.0743802,0.0413223];
vector25=[0.0184911,0.035503,0.0510355,0.0650888,0.0776627,0.0887574,0.0983728,0.106509,0.113166,0.118343,0.122041,0.12426,0.125,0.12426,0.122041,0.118343,0.113166,0.106509,0.0983728,0.0887574,0.0776627,0.0650888,0.0510355,0.035503,0.0184911];
vector100=[0.00490148,0.00970493,0.0144104,0.0190177,0.0235271,0.0279384,0.0322517,0.036467,0.0405843,0.0446035,0.0485247,0.0523478,0.0560729,0.0597,0.0632291,0.0666601,0.0699931,0.0732281,0.0763651,0.079404,0.0823449,0.0851877,0.0879326,0.0905794,0.0931281,0.0955789,0.0979316,0.100186,0.102343,0.104402,0.106362,0.108225,0.109989,0.111656,0.113224,0.114695,0.116067,0.117341,0.118518,0.119596,0.120576,0.121459,0.122243,0.122929,0.123517,0.124007,0.1244,0.124694,0.12489,0.124988,0.124988,0.12489,0.124694,0.1244,0.124007,0.123517,0.122929,0.122243,0.121459,0.120576,0.119596,0.118518,0.117341,0.116067,0.114695,0.113224,0.111656,0.109989,0.108225,0.106362,0.104402,0.102343,0.100186,0.0979316,0.0955789,0.0931281,0.0905794,0.0879326,0.0851877,0.0823449,0.079404,0.0763651,0.0732281,0.0699931,0.0666601,0.0632291,0.0597,0.0560729,0.0523478,0.0485247,0.0446035,0.0405843,0.036467,0.0322517,0.0279384,0.0235271,0.0190177,0.0144104,0.00970493,0.00490148];

xvec10=zeros(10,1);
xvec25=zeros(25,1);
xvec100=zeros(100,1);

for i=1:10
    xvec10(i)=(i+1)/11;
end


for i=1:25
    xvec25(i)=(i+1)/26;
end

for i=1:100
    xvec100(i)=(i+1)/101;
end

plot(xvec10,vector10,'r*');
xlabel('(i+1)/(n+1)','FontSize', 14);
ylabel('Final vector xi','FontSize', 14);
title('N=10','FontSize', 14);
figure
plot(xvec25,vector25,'r*');
xlabel('(i+1)/(n+1)','FontSize', 14);
ylabel('Final vector xi','FontSize', 14);
title('N=25','FontSize', 14);
figure
plot(xvec100,vector100,'r*');
xlabel('(i+1)/(n+1)','FontSize', 14);
ylabel('Final vector xi','FontSize', 14);
title('N=100','FontSize', 14);