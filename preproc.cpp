

void read_setup_file(const std::string &filename, freestream &free_stream, material &mat, boundary_marker_list &bound, cellist_2D &clist, std::string &meshfile, simparam &sim)
{

    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error: could not open " << filename << "\n";
        std::exit(1);
    }

    file >> meshfile;
    file >> free_stream.Mach;
    file >> free_stream.AoA;
    // file>>free_stream.rho_inf;
    file >> free_stream.p_inf;
    file >> free_stream.T_inf;
    file >> mat.gamma;
    file >> mat.R;
    file >> sim.CFL;
    file >> sim.tol;
    file >> sim.MAX_ITER;
    file >> sim.print_interval;
    file >> sim.write_interval;
    file >> sim.output_filename;

    std::cout << "|| Reading grid..... ||" << "\n";

    read_grid(clist, bound, meshfile);


    // loop over boundary markers
    for (int i = 0; i < bound.nmark; i++)
    {
        std::string markername;
        int temp_type=0;

        file>>markername>>temp_type; 

       // std::cout<<markername<<","<<temp_type<<std::endl;

        if(markername==bound.marker_list[i].marker_tag)
        {
            bound.marker_list[i].btype = temp_type;
        }
        
    }

    file.close();
}



