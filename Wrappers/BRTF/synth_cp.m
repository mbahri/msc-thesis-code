U1 = rand(120,20); U2 = rand(110,20); U3 = rand(100,20);
O = double(ktensor(ones(20,1), U1, U2, U3));

W = randi([0,1], 120, 110, 100);
E = sign(rand(120, 110, 100) - 0.5) .* W;

X = O + E;