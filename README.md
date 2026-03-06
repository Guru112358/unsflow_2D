This is a simple cell centered Finite Volume Euler solver on triangular meshes using the SU2 mesh format(filename.su2)  as an input .
It uses Rusanov Flux with a MUSCL reconstruction and offers a choice of cell based Barth and Venkatakrishnan limiters, it has been validated here using the RAE 2822 transonic airfoil case with the Benchmark being an SU2 solution on the same grid .

The solver currently onyl supports Euler and far field boundary conditions and can be run after compilation simply by the command " ./unsflow setup.inp" where setup.inp is the setup file that has all the relevant information on how to go about setting parameters.

THe plots are all handled by gnuplot and the post processing is handled by paraview (output of the program is in vtk format) .


<img width="1920" height="947" alt="0 725_2 31" src="https://github.com/user-attachments/assets/f08240a2-7c0d-4305-bf16-ecfccb10e1c4" />
