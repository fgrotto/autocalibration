# Simple Autocalibration Procedure

Assumptions: all the images have the same internal parameters K

 Steps to perform:
- [ ] Prepare the data structure (nodes for pair of images i,j with 3D point and (u,v) in the images called corrispondence points). Use the script pose_driver.m provided 
- [ ] Then we can compute F (fundamental matrix) with 8 points algorithm (implemented in the toolbox) by using points from left and right image in both direction
- [ ] Use the script provided in the autocal internal_driver to get an estimation of K. We call it K0
- [ ] Use autocalibration procedure medonca-cipolla to start and get the estimated K
- [ ] From the estimated K compare it with the real K (the one used in pose_driver.m to map 3d to the image pixels (u,v)
- [ ] From the estimated K compute R|t relative to each pair of image (i,j). So we will have the rotation R and translation t of j wrt of i.
- [ ] By triangulation using ppm (with K,R,t estimated) project the 2d points to 3d and get a cloud of points.
- [ ] Register the cloud of points with ICP to compare it with the original (basically the problem is that all the estimated R|t are of j wrt of i not the same reference frame as the original cloud point but we can move everthing to i or by concatenation)

Suggestions:
- Compute F fundamental matrix on images with sufficient overlap (take 10 points).
- Autocalibration uses the fundamental matrixes computed, see if we need to get all of them.
- Maybe consider outliners with RANSAC (the prof didn't mention it but maybe we can check it).
- Maybe consider other autocalibration procedures.
