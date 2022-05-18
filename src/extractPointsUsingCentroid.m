function close_points = extractPointsUsingCentroid(points, distance)
   m = mean(points);
   
   index = 1;
   done = false;
   for i = 1:size(points,1)
        if norm(points(i,:)-m) < distance
            close_points(index,:) = points(i,:);
            index = index+1;
            done = true;
        end
   end
   
   if ~done
       close_points = [];
   end
end