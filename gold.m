function code = gold(code_serial) 
% convention is {31, 63, 127, 511, 1023, 2047} correspond to {1,2,3,4,5,6}
codes = [31, 63, 127, 511, 1023, 2047];
M = codes(code_serial);
preferred_pairs_1 = [5 2 0, 
6 1 0, 
7 3 0, 
9 4 0, 
10 3 0, 
11 2 0];
preferred_pairs_2 = [5 4 3 2 0,
6,5,2,1,0,
7,3,2,1,0,
9,6,4,3,0,
10,8,3,2,0,
11,8,5,2,0];

pn = comm.PNSequence('Polynomial', preferred_pairs_1(code_serial, :), 'SamplesPerFrame', 2*codes(code_serial), 'InitialConditions', [0 0 0 0 1]);
x1 = step(pn);

x1 = x1(1:M);

pn = comm.PNSequence('Polynomial', preferred_pairs_2(code_serial, :), 'SamplesPerFrame', 2*codes(code_serial), 'InitialConditions', [0 0 0 0 1]);
x2 = step(pn);

%[x2(1:31) x2(32:62)]
x2 = x2(1:M);

% x1 remaining fixed, x2 is shifting crcularly.
codes = zeros(M+2, M);
codes(1,:) = x1;
codes(2,:) = x2;
for i=3:M+2
    codes(i,:) = xor(circshift(x2,i-2), x1);
end

correlation = zeros(M+2, M+2);

codes = 2 * (codes>0) - 1;

for i=1:M+2
  for j=1:M+2
    correlation(i,j) = max(abs(xcorr(codes(i,:),codes(j,:))));
  end
end

min_cross_corr = ones(1,M+2);
min_id = ones(1,M+2);

for i=1:M+2
  [min_cross_corr(i), min_id(i)] = min(correlation(i,:));
end

[k, k_idx] = min(min_cross_corr);

stem(xcorr(codes(k_idx,:), codes(min_id(k_idx),:)));
fprintf('%d %d\n', k_idx, min_id(k_idx));

% AutoCorr, CrossCorr: use xcorr() : nC2 : get min(max) choose from this, and put.
end
