In this project, I adapted the tensor voting framework originally written by Trevor Linton in Matlab, into Java, to be utilized in a cybersecurity app.

Tensor voting uses the Gesalt principles of proximity and continuation to infer salient structures from a set of points. 
Tensor voting is a combination of tensor calculus, for point representation, and non-linear voting to transfer data. 
All points are represented by a tensor and a saliency value. 
Saliency tells the "confidence" of a point, or the strength a point will have on influencing its neighbors. 
Each coordinate point is composed of two seperate tensors, a ball and a stick tensor. 
A stick tensor is a 2x2 elipsoid that represents the dirrection the point is aimed, and its strength. 
A ball tensor contributes to the strength of a point when performing non-linear voting, but does not provide a specific direction.

After tensor values are calulated, non-linear voting is performed in which the tensor of each point is impacted by the tensors of every other point in its defined neighborhood. 
This causes points in close proximity to each other to become grouped, creating lines between them. 
The only paramater in tensor voting that can be manually changed is the scale, which decides the radius of neighborhoods used in non-linear voting. 
As the scale is increased, the framework becomes less responsive to noise, but gives us detail. 