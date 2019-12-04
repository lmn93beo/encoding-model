function spline_basis = make_delta_splines(T, filename)
spline_basis = eye(T);
save(fullfile('Splines', filename), 'spline_basis');

