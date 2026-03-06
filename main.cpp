#include<iostream>
#include<fstream>
#include<array>
#include<cmath>
#include<vector>
#include<string>
#include <sstream>
#include "materials.h"

#include "flowvar.cpp"
#include "simparam.h"
#include "geometry.h"
#include "grid.cpp"
#include "fvclass.h"
#include "euler_solve.cpp"
#include "postproc.h"
#include "init.cpp"
#include "preproc.cpp"


int main(int argc, char* argv[])
{

    material mat;
    simparam sim;

    finite_volume_2D fv;

    freestream free_stream;

    std::string meshfile;

    read_setup_file(argv[1],free_stream,mat,fv.boundary,fv.clist,meshfile,sim);


    std::cout<<"||  computing  geometry parameters ..... ||"<<"\n";


    compute_centroid_and_grid_volume_2D(fv.clist);

    std::cout<<"||  Forming faces  ..... ||"<<"\n";

    form_faces(fv.clist,fv.F);

    std::cout<<"||  Finding neighbours and setting boundary faces  ..... ||"<<"\n";


    find_neighbours(fv.clist,fv.F);

    set_boundary_faces(fv.boundary,fv.F);

    free_stream= compute_freestream(free_stream,mat); 
 

    init_freestream(fv.clist,free_stream,mat);


    std::cout<<"||  starting UNSFLOW: - ||"<<"\n";


        compute_residual(fv.clist,fv.boundary,fv.F, mat,sim);  

        apply_boundary_conditions(fv.boundary,fv.clist , fv.F, free_stream, mat);

        auto res0 = compute_L2_residual(fv.clist);

        bool loop_switch=true;

        long int loop_counter=0;

        std::array<double, 4> res;

        std::fstream res_record;

        res_record.open("convergence_history.csv",std::ios::out);

    while(loop_switch)
    {

        compute_residual(fv.clist,fv.boundary,fv.F, mat,sim);  

        apply_boundary_conditions(fv.boundary,fv.clist , fv.F, free_stream, mat);

        solution_update(fv.clist);


          
        for (size_t i = 0; i < fv.clist.cell_list.size(); ++i)
        {
            compute_primitives(fv.clist.cell_list[i].prim,fv.clist.cell_list[i].cons,mat);
        }

       	
        res =compute_L2_residual(fv.clist);

        if(loop_counter%sim.print_interval==0)
        {
            compute_pressure_coefficient("cp.csv",fv.boundary,fv.F,fv.clist,mat,free_stream,1);

            std::cout << "Iter " << loop_counter << " | R = ["
          << res[0]/res0[0] << ", "
          << res[1]/res0[1] << ", "
          << res[2]/res0[2] << ", "
          << res[3]/res0[3] << "]"
         // <<"  | Cl= " << fv.boundary.marker_list[1].coeffs[1]   (they seem to be a bit broken at the moment so would not suggest trusting the values)
          //<<", Cd= " << fv.boundary.marker_list[1].coeffs[0] <<" |"
          << std::endl;

          res_record<<loop_counter << ","
          << res[0]/res0[0] << ", "
          << res[1]/res0[1] << ", "
          << res[2]/res0[2] << ", "
          << res[3]/res0[3] << ", "
          //<< fv.boundary.marker_list[1].coeffs[1]
         // <<", " << fv.boundary.marker_list[1].coeffs[0]
			  <<std::endl;
             //write_vtk_2d("output.vtk",fv.clist.cell_list);

        
        }

        if(loop_counter%sim.write_interval==0)
        {
                write_vtk_2d(sim.output_filename,fv.clist.cell_list);
                compute_pressure_coefficient("cp.csv",fv.boundary,fv.F,fv.clist,mat,free_stream,1);
        }

        
        if((res[0]/res0[0]<sim.tol&&(loop_counter!=0)))
        {
        std::cout<<"converged to a tolerance of: "<<sim.tol<<" in "	<<loop_counter<<" Iterations"<<"\n";
        loop_switch=false;
        }

        if(loop_counter==sim.MAX_ITER)
        {
            std::cout<<"Exiting, Maximum iterations reached: "<<sim.MAX_ITER<<"\n";
	        loop_switch=false;
        }

   // std::cout<<i<<std::endl;
         loop_counter++;

    }

    write_vtk_2d(sim.output_filename,fv.clist.cell_list);

    compute_pressure_coefficient("cp.csv",fv.boundary,fv.F,fv.clist,mat,free_stream,1);

    std::cout<<"| Cl= " << fv.boundary.marker_list[1].coeffs[1]<<", Cd= " << fv.boundary.marker_list[1].coeffs[0]
    <<" |"<<std::endl;


    return 0;
 
}
