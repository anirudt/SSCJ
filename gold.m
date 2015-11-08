clear all
close all

pn = comm.PNSequence('Polynomial', [5 4 3 1 0], 'SamplesPerFrame', 62, 'InitialConditions', [0 0 0 0 1]);
x1 = step(pn);

%[x1(1:31) x1(32:62)]
x1 = x1(1:31);

pn = comm.PNSequence('Polynomial', [5 4 3 2 0], 'SamplesPerFrame', 62, 'InitialConditions', [0 0 0 0 1]);
x2 = step(pn);

%[x2(1:31) x2(32:62)]
x2 = x2(1:31);

% x1 remaining fixed, x2 is shifting crcularly.
codes = zeros(33, 31);
codes(1,:) = x1;
codes(2,:) = x2;
for i=3:33
    codes(i,:) = xor(circshift(x2,i-2), x1);
end

correlation = zeros(33, 33);

codes = 2 * (codes>0) - 1;

for i=1:33
  for j=1:33
    correlation(i,j) = max(abs(xcorr(codes(i,:),codes(j,:))));
  end
end

min_cross_corr = ones(1,33);
min_id = ones(1,33);

for i=1:33
  [min_cross_corr(i), min_id(i)] = min(correlation(i,:));
end

[k, k_idx] = min(min_cross_corr);

stem(xcorr(codes(k_idx,:), codes(min_id(k_idx),:)));
fprintf('%d %d\n', k_idx, min_id(k_idx));

% AutoCorr, CrossCorr: use xcorr() : nC2 : get min(max) choose from this, and put.
