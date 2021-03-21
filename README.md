# Design Optimization of a Pneumatic Soft Robotic Actuator 
Soft robotic devices offer superior solutions for biomechanical engineering problems at the human-machine interface due to their compliant nature. To date limited research has explored the topic of surface manipulation which can be helpful in the prevention of pressure injuries. The codes in this repository are concerned with the shape optimzation of air chambers to achive maximal horizontal motion with minimal vertical motion from both an analytical model and a finite element based approach. We explore both the firefly algorithm and the deep determinisitc policy gradient (DDPG) coupled with a dervied analytical model, and DDPG couple with a finite element program in this repository.





The model based approach considers the shown actuator as a system of cantilever beams connected by spring elements. The governing eqations for this as given in our paper and are used as the cost function in the firefly approach. In the DDPG based approach a vanilla gym structure is used with the model embedded into the step function. In the finite element based approach the system is modeled as a quadrilateral whose verticies can exist within fixed bounds in the semi-unsupervised learning problem. We reward larger horizontal motion and penalize vertical motion in this approach.

# Requirements
The DDPG program requires pyansys and Stable Baselines for use. Please see the Gym instructions text file for an explaination of where to place the environment files.

# Publication
Our publication for this work is still under review and will be updated soon!
