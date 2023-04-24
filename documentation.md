## Abstract

Many novel biomaterials are recently investigated for use in spinal fusion surgery, especially in lumbar interbody fusion. The X-ray microCT as a tool is widely used for evaluating how successfully those biomaterials can perform a vertebral fusion. 
  
  Our methodology presents an automatic approach for pig’s spinal fusion evaluation in 3D. The proposed approach is based on the determination of the vertebral fused area, which reflects the fusion quality.The calculation is based on the detection of the fused area and area of facies intervertebralis, so the percentual representation of the vertebral joint can be determined.
  

![2C16 - Copy2](https://user-images.githubusercontent.com/41157503/234028034-993c98f4-bdf2-45ee-a958-9bded2b10dbc.png)

Figure: X-ray micro CT measurement of two fused pig vertebrae. Red: Vertebral body area (expanded and depicted in red for visualization). Green: Detected area of vertebral fusion.


## Code functions:
The algorithm is located in the folder "Main"

**1. Main.m** Main function of the software. Inputs and outputs of the algorithm are inserted here.

**2. detect_concave_and_convex_points.m** Detection of convex and concave points in order to detect the location of both fused vertebrae in the dataset.

**3. convex_or_concave.m** Returns the convex and concave points in vertebrae boundaries - decision if the point is convex or concave

**4. find_quadrilaterals.m** Fitting Quadrilaterals to approximate upper and bottom vertebrae.

**5. compute_rotated_ellipses.m** The function returns coordinates of ellipses aproximating upper and bottom vertebrae

**6. intersection_of_ellipses.m** Returns coordinates of intersections of two ellipses approximating both vertebrae.

**7. intersection_of_lines.m** Returns coordinates of intersection of two lines defined by centroids of ellipses and intersections of ellipses. 

**8. find_roi.m** Returns the area of ROI, according to the coordinates of ellipses and lines intersection. In this ROI is located the fused area.

**9. apply_watershed.m** Function which separates both vertebrae in order to find the vertebral body and fusion areas.

**10. -area_ratio.m** Final computation of the ratio between vertebral bodies surfaces and fused area in order to evaluate, how sucessfully both vertebrae fused.

**-mia_curve_tangent.m** Supportive function to detect_concave_and_convex_points.m. Function to extract curves from binary edge map, if the endpoint of a contour is nearly connected to another endpoint, fill the gap and continue the extraction.

**-mia_bresenham.m** Supportive function to convex_or_concave.m. Function creates a line beetween two point based on bresenham  algorithm.

** More description is provided in individual code files.

## Testing dataset:
Available here: Data/Link_to_repository.url

**Email me: Jakub.Laznovsky@ceitec.vutbr.cz**
