function error = cost_kruppas_method(Fs, k)
    % Compute kruppas method by passing Fs the fundamental matrixes and x
    % and initial estimation of the intrinsic parameters (in vector form)
    % to allow the calculation to converge to the correct value.
    K=[k(1) k(2) k(3); 0 k(4) k(5); 0 0 1];

    w = K * K';
    error = [];

    for i = 1:size(Fs,1)
        for j = i+1:size(Fs,2)
            if ~isempty(Fs{i,j})
                Kleft = Fs{i,j} * w * Fs{i,j}';
                Kleft = Kleft/norm(Kleft,'fro');

                Kright = epipole(Fs{i,j}) * w * epipole(Fs{i,j})';
                Kright = Kright/norm(Kright, 'fro');

                % compute the kruppas error
                Kdiff = Kleft - Kright;
                error = [error Kdiff(1,:) Kdiff(2,2:3)];
            end
        end
    end
end

