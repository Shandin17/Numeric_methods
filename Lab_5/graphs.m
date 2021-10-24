T = readmatrix('Characteristics.dat');
X = readmatrix('X.dat');

for i = 1: length(X(:,1))
plot(X(i,:),T(i,:));
hold on
end