function pts = compute_points_from_line(lines, imageSize)
    % Adapted from matlab code
    nPts = size(lines, 1);
    pts = coder.nullcopy(-ones(nPts, 4));
    firstRow = 0.5;
    firstCol = 0.5;

    lastRow = firstRow + imageSize(1);
    lastCol = firstCol + imageSize(2);

    % Loop through all lines and compute the intersection points of the lines
    % and the image border.
    for iLine = 1:nPts
      a = lines(iLine, 2);
      b = lines(iLine, 1);
      c = lines(iLine, 3);

      endPoints = zeros([4, 1]);
      iPoint = ones(1);
      % Check for the intersections with the left and right image borders
      % unless the line is vertical.
      if abs(a) > 0.001

        % Compute and check the intersection of the line and the left image
        % border. 
        row = - (b * firstCol + c) / a;
        if row>=firstRow && row<=lastRow
          endPoints(iPoint:iPoint+1) = [row; firstCol];
          iPoint = iPoint + 2;
        end

        % Compute and check the intersection of the line and the right image
        % border. 
        row = - (b * lastCol + c) / a;
        if row>=firstRow && row<=lastRow
          endPoints(iPoint:iPoint+1) = [row; lastCol];
          iPoint = iPoint + 2;
        end
      end

      % Check for the intersections with the top and bottom image borders
      % unless the line is horizontal.
      if abs(b) > 0.001
        % If we have not found two intersection points, compute and check the
        % intersection of the line and the top image border. 
        if iPoint < 4
          col = - (a * firstRow + c) / b;
          if col>=firstCol && col<=lastCol
            endPoints(iPoint:iPoint+1) = [firstRow; col];
            iPoint = iPoint + 2;
          end
        end

        % If we have not found two intersection points, compute and check the
        % intersection of the line and the bottom image border. 
        if iPoint < 4
          col = - (a * lastRow + c) / b;
          if col>=firstCol && col<=lastCol
            endPoints(iPoint:iPoint+1) = [lastRow; col];
            iPoint = iPoint + 2;
          end
        end
      end

      % If the line does not intersect with the image border, set the
      % intersection to -1; 
      for i = iPoint: 4
        endPoints(i) = -1;
      end

      pts(iLine, :) = endPoints([2,1,4,3]);
    end
end
