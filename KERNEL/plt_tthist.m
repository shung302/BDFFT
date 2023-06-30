figure
p=load('pttvector.dat');
s=load('sttvector.dat');
h1 = histogram(p,16);
%h1.Normalization = 'probability';
h1.BinWidth = 0.25;
hold on
h2=histogram(s,16);
%h2.Normalization = 'probability';
h2.BinWidth = 0.25;
h2.FaceColor = [0 0.5 0.5];
h2.EdgeColor = 'r';
axis([-0.8 0.8 0 100])
xlabel('P and S Residual (s)');
ylabel('Count');
