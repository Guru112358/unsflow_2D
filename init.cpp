   
   
   void init_freestream(cellist_2D &clist,freestream &free_stream,material &mat)
   {

   for (auto &cell : clist.cell_list)
    {

        cell.cons.U[0] = free_stream.rho_inf;
        cell.cons.U[1] = free_stream.rho_inf * free_stream.u_inf;
        cell.cons.U[2] = free_stream.rho_inf * free_stream.v_inf;
        cell.cons.U[3] = free_stream.rho_inf * free_stream.E_inf;

        cell.prim.rho = cell.cons.U[0];
        cell.prim.u   = cell.cons.U[1] / cell.prim.rho;
        cell.prim.v   = cell.cons.U[2] / cell.prim.rho;

        double kinetic = 0.5 * (cell.prim.u*cell.prim.u +
                                cell.prim.v*cell.prim.v);

        double E = cell.cons.U[3] / cell.prim.rho;

        cell.prim.p = (mat.gamma - 1.0)
                      * cell.prim.rho
                      * (E - kinetic);

    }
}