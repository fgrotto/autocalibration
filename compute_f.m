function [Fs, S] = compute_f(S, directory, indexes, num_points)
    ind = 1;
    for x = 1:length(indexes)
        for y = 1:length(indexes)
           i = indexes(x);
           j = indexes(y);

           if i == j 
               continue
           end
           I_left = [directory S{i,j}.name_view_i];
           I_right = [directory S{i,j}.name_view_j];
            
           G = [S{i,j}.R_i S{i,j}.t_i; 0 0 0 1];
           P1 = S{i,j}.K_i * [1 0 0 0
                    0 1 0 0
                    0 0 1 0] * G;
                
           G = [S{i,j}.R_j S{i,j}.t_j; 0 0 0 1];
           P2 = S{i,j}.K_j * [1 0 0 0
                    0 1 0 0
                    0 0 1 0] * G;
                
           % Simple test with triangulation to check the result
           % X = triang(P1,P2, [S{i,j}.uv_i(1:1) S{i,j}.vv_i(1:1)]', [S{i,j}.uv_j(1:1) S{i,j}.vv_j(1:1)]')
           % S{1,2}.points(1,:)

           if (size(S{i,j}.uv_i,1) > num_points) && (size(S{i,j}.uv_j,1) > num_points)
               leftP = [S{i,j}.uv_i(:) S{i,j}.vv_i(:)]';
               rightP = [S{i,j}.uv_j(:) S{i,j}.vv_j(:)]';

               [F_expected, ~] = fund(P1,P2);
               [F, left_in, right_in, inliers] = fundamental_matrix(I_left, I_right, leftP, rightP, num_points);
               S{i,j}.inliers = inliers;

               fprintf('Fundamental nonlin Smps error views (%0.2g, %0.2g):\t %0.5g \n', i,j, rmse(sampson_fund(F,left_in,right_in)));
               if (rmse(sampson_fund(F,leftP,rightP)) < 1)
                   Fs(:,:,1,ind) = F;
                   ind=ind+1;
               end
           end
        end
    end
end