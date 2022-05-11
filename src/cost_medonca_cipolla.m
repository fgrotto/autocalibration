function C = cost_medonca_cipolla(F, X)
    %   Inputs
    %       F - Fundamenta matrixes given images
    %       X - Initial estimation of instrinsic camera parameters (vector
    %       form)
    %   Output
    %       E - Computed Cost
    
    % Get initial K estimation
    K = [X(1) X(2) X(3); 0 X(4) X(5); 0 0 1];

    % Initialize Cost
    C = [];
    
    % Obtain the number of images
    n = size(F,1)*size(F,2);
    
    % By matching points pairwise it is possible to find
    % n(n-1)/2 fundamental matrixes
    Den = n*(n-1)/2;

    for i = 1:size(F,1)
        for j = (i+1):size(F,2)
            if ~isempty(F{i,j})
                % Essential matrix from fundamental matrix
                E = K' * F{i,j} * K;

                % SVD of essential matrix
                [~,D,~] = svd(E);

                % Non zero singular values of 
                sigma1 = D(1,1);
                sigma2 = D(2,2);

                C1 = 1/Den * (sigma1 - sigma2)/sigma2;
                C = [C C1];
            end
        end
    end
end
