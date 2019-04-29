function [v, H, R] = associateMeasurements(x, P, Z, R, M, g)
% [v, H, R] = associateMeasurements(x, P, Z, R, M, g) returns a set of
% innovation vectors and associated jacobians and measurement covariances
% by matching line features by Mahalanobis distance.


nMeasurements = size(Z, 2);
nMapEntries = size(M, 2);

d = zeros(nMeasurements, nMapEntries);
v = zeros(2, nMeasurements * nMapEntries);
H = zeros(2, 3, nMeasurements * nMapEntries);

for i = 1:nMeasurements
    for j = 1:nMapEntries
        [z_priori, H(:, :, j+(i-1)*nMapEntries)] = measurementFunction(x, M(:, j));
        v(:, j+(i-1)*nMapEntries) = Z(:, i) - z_priori;
        W = H(:, :, j+(i-1)*nMapEntries) * P * H(:, :, j+(i-1)*nMapEntries)' + R(:, :, i);
        d(i, j) = v(:, j + (i-1)*nMapEntries)' * inv(W) * v(:, j+(i-1)*nMapEntries);
    end
end

[minima, map_index] = min(d');
[measurement_index] = find(minima < g^2);
map_index = map_index(measurement_index);

v = v(:, map_index + (measurement_index-1)*nMapEntries);
H = H(:, :, map_index+(measurement_index-1)*nMapEntries);
R = R(:, :, measurement_index);