void compute_residual(cellist_2D &C, boundary_marker_list &boundary, facelist_2D &F, material &mat, simparam &sim);
void solution_update(boundary_marker_list &boundary, facelist_2D &F);
void apply_boundary_conditions(boundary_marker_list &boundary, facelist_2D &F);
void compute_Lax_Friedrichs_flux(face2D &current_face, cellist_2D &C, flux &result, material &mat);
void compute_fluxes_cell_center(flux_2D &flux, primitive &var, material &mat);
void compute_normal_flux(face2D &current_face, flux &Fn, primitive &var, material &mat);
double compute_cell_time_step(cell2D &C, facelist_2D &F, double CFL, material &mat);
void compute_Lax_Friedrichs_flux_freestream(face2D &current_face, cellist_2D &C, cell2D &ghost_neighbour, flux &result, material &mat, freestream &free_stream);
double compute_lambda(double un_L, double un_R, double P_L, double P_R, double rho_L, double rho_R, double gamma);
void reset_residuals(cellist_2D &C);
void apply_freestream(boundary_marker_list &boundary, cellist_2D &C, facelist_2D &F, freestream &free_stream, material &mat, int marker_index);
void apply_Euler_wall(boundary_marker_list &boundary, cellist_2D &C, facelist_2D &F, freestream &free_stream, material &mat,int marker_index);
std::array<double, 4> compute_L2_residual(cellist_2D &C);
void compute_Gauss_Green_gradient(int cell_index ,facelist_2D &F, cellist_2D &C);
void compute_conserved_cell_center(primitive &var, conserved &cons, material &mat);
void correct_Gauss_Green_gradient(int cell_index ,facelist_2D &F,cellist_2D &C);
void compute_barth_limiter_cell(int cell_index, facelist_2D &F, cellist_2D &C);
primitive reconstruct_face_primitive_limited(const cell2D &C,const point2D &xf);





inline double barth_limiter(double dq, double qC, double qmin, double qmax)
    {
        double phi = 1.0;

        if (dq > 0.0)
            phi = std::min(1.0, (qmax - qC) / dq);
        else if (dq < 0.0)
            phi = std::min(1.0, (qmin - qC) / dq);

        return std::max(0.0, phi);
    };


inline double venkat_limiter(double dq,double qC,double qmin,double qmax, double eps2)
{
    if (dq == 0.0) return 1.0;

    double delta;

    if (dq > 0.0)
        delta = qmax - qC;
    else
        delta = qmin - qC;

    double num = delta*delta + 2.0*delta*dq + eps2;
    double den = delta*delta + 2.0*dq*dq + delta*dq + eps2;

    if (den == 0.0) return 1.0;

    double phi = num / den;

    if (phi < 0.0) phi = 0.0;
    if (phi > 1.0) phi = 1.0;

    return phi;
}

void compute_barth_limiter_cell(int cell_index, facelist_2D &F, cellist_2D &C)
{
    cell2D &cell = C.cell_list[cell_index];

    double phi_r = 1.0, phi_u = 1.0, phi_v = 1.0, phi_p = 1.0;

    double rC = cell.prim.rho;
    double uC = cell.prim.u;
    double vC = cell.prim.v;
    double pC = cell.prim.p;

    double rmin = rC, rmax = rC;
    double umin = uC, umax = uC;
    double vmin = vC, vmax = vC;
    double pmin = pC, pmax = pC;

    for (int f = 0; f < 3; f++)
    {
        int nb = cell.neighbour_cell_number[f];
        if (nb < 0) continue;

        const cell2D &N = C.cell_list[nb];

        rmin = std::min(rmin, N.prim.rho);
        rmax = std::max(rmax, N.prim.rho);

        umin = std::min(umin, N.prim.u);
        umax = std::max(umax, N.prim.u);

        vmin = std::min(vmin, N.prim.v);
        vmax = std::max(vmax, N.prim.v);

        pmin = std::min(pmin, N.prim.p);
        pmax = std::max(pmax, N.prim.p);
    }

    for (int f = 0; f < 3; f++)
    {
        int fid = cell.face_id_indices[f];
        const face2D &face = F.face_list[fid];

        const point2D &xf = face.face_centroid;

        double dx = xf.pos[0] - cell.centroid.pos[0];
        double dy = xf.pos[1] - cell.centroid.pos[1];

        double dr = cell.Grad.rho[0]*dx + cell.Grad.rho[1]*dy;
        double du = cell.Grad.u[0]*dx   + cell.Grad.u[1]*dy;
        double dv = cell.Grad.v[0]*dx   + cell.Grad.v[1]*dy;
        double dp = cell.Grad.p[0]*dx   + cell.Grad.p[1]*dy;

        phi_r = std::min(phi_r, barth_limiter(dr, rC, rmin, rmax));
        phi_u = std::min(phi_u, barth_limiter(du, uC, umin, umax));
        phi_v = std::min(phi_v, barth_limiter(dv, vC, vmin, vmax));
        phi_p = std::min(phi_p, barth_limiter(dp, pC, pmin, pmax));
    }

    cell.limiter_rho = phi_r;
    cell.limiter_u   = phi_u;
    cell.limiter_v   = phi_v;
    cell.limiter_p   = phi_p;
}

void compute_venkat_limiter_cell(int cell_index,facelist_2D &F,cellist_2D &C)
{
    cell2D &cell = C.cell_list[cell_index];

    double phi_r = 1.0, phi_u = 1.0, phi_v = 1.0, phi_p = 1.0;

    double rC = cell.prim.rho;
    double uC = cell.prim.u;
    double vC = cell.prim.v;
    double pC = cell.prim.p;

    double rmin = rC, rmax = rC;
    double umin = uC, umax = uC;
    double vmin = vC, vmax = vC;
    double pmin = pC, pmax = pC;

    for (int f = 0; f < 3; f++)
    {
        int nb = cell.neighbour_cell_number[f];
        if (nb < 0) continue;

        const cell2D &N = C.cell_list[nb];

        rmin = std::min(rmin, N.prim.rho);
        rmax = std::max(rmax, N.prim.rho);

        umin = std::min(umin, N.prim.u);
        umax = std::max(umax, N.prim.u);

        vmin = std::min(vmin, N.prim.v);
        vmax = std::max(vmax, N.prim.v);

        pmin = std::min(pmin, N.prim.p);
        pmax = std::max(pmax, N.prim.p);
    }

    double h = 1 ; //following SU2 documentation here, K=0
    double K = 0.01;                   
    double eps2 = std::pow(K*h, 3.0);

    for (int f = 0; f < 3; f++)
    {
        int fid = cell.face_id_indices[f];
        const face2D &face = F.face_list[fid];

        double dx = face.face_centroid.pos[0] - cell.centroid.pos[0];
        double dy = face.face_centroid.pos[1] - cell.centroid.pos[1];

        double dr = cell.Grad.rho[0]*dx + cell.Grad.rho[1]*dy;
        double du = cell.Grad.u[0]*dx   + cell.Grad.u[1]*dy;
        double dv = cell.Grad.v[0]*dx   + cell.Grad.v[1]*dy;
        double dp = cell.Grad.p[0]*dx   + cell.Grad.p[1]*dy;

        phi_r = std::min(phi_r, venkat_limiter(dr, rC, rmin, rmax, eps2));
        phi_u = std::min(phi_u, venkat_limiter(du, uC, umin, umax, eps2));
        phi_v = std::min(phi_v, venkat_limiter(dv, vC, vmin, vmax, eps2));
        phi_p = std::min(phi_p, venkat_limiter(dp, pC, pmin, pmax, eps2));
    }

    cell.limiter_rho = phi_r;
    cell.limiter_u   = phi_u;
    cell.limiter_v   = phi_v;
    cell.limiter_p   = phi_p;
}

void compute_Gauss_Green_gradient(int cell_index ,facelist_2D &F,cellist_2D &C)
{

    cell2D &cell= C.cell_list[cell_index];
    
 
    for(int f=0;f<3;f++)
    {

        int fid =cell.face_id_indices[f];

        face2D &current_face = F.face_list[fid];
        
        primitive phi_face_temp = cell.prim;

        int nb = cell.neighbour_cell_number[f];

        if(nb>=0)// not boundary cell
        {

        cell2D  &neighbour = C.cell_list[nb];

        double d_owner_neighbour = std::sqrt
        (
           std::pow( neighbour.centroid.pos[0] - cell.centroid.pos[0], 2) +
           std::pow( neighbour.centroid.pos[1] - cell.centroid.pos[1], 2)
        );
        
        double d_face_neighbour = std::sqrt
        (
           std::pow( current_face.face_centroid.pos[0] - neighbour.centroid.pos[0], 2) +
           std::pow( current_face.face_centroid.pos[1] - neighbour.centroid.pos[1], 2)
        );
        

        //refer :CFD textbook by Moukalled and Darwish
         double gc = d_face_neighbour/d_owner_neighbour;

            phi_face_temp.rho =  gc*(cell.prim.rho ) + (1.0-gc)*(neighbour.prim.rho);
            phi_face_temp.u =    gc*(cell.prim.u) + (1.0-gc)*(neighbour.prim.u);
            phi_face_temp.v =    gc*(cell.prim.v) + (1.0-gc)*(neighbour.prim.v);
            phi_face_temp.p =    gc*(cell.prim.p) + (1.0-gc)*(neighbour.prim.p);
           
        }

         double factor = current_face.len / cell.cell_area;     
         double nx = current_face.n.pos[0];
         double ny = current_face.n.pos[1];

         //checking normals for sign consistency
         double sign = (current_face.owner == cell_index) ? 1.0 : -1.0;
         nx *= sign;
         ny *= sign;

         cell.Grad.rho[0] +=  phi_face_temp.rho* nx * factor ;
         cell.Grad.rho[1] +=  phi_face_temp.rho* ny * factor ;

         cell.Grad.u[0] +=  phi_face_temp.u* nx * factor ;
         cell.Grad.u[1] +=  phi_face_temp.u* ny * factor ;

         cell.Grad.v[0] +=  phi_face_temp.v* nx * factor ;
         cell.Grad.v[1] +=  phi_face_temp.v* ny * factor ;

         cell.Grad.p[0] +=  phi_face_temp.p* nx * factor ;
         cell.Grad.p[1] +=  phi_face_temp.p* ny * factor ;

    }

   //  std::cout<< cell.Grad.u[0]<<","<< cell.Grad.u[1]<<std::endl;

}


primitive reconstruct_face_primitive_limited(const cell2D &C,const point2D &xf)
{
    primitive W = C.prim;

    double dx = xf.pos[0] - C.centroid.pos[0];
    double dy = xf.pos[1] - C.centroid.pos[1];

    double dr = C.Grad.rho[0]*dx + C.Grad.rho[1]*dy;
    double du = C.Grad.u[0]  *dx + C.Grad.u[1]  *dy;
    double dv = C.Grad.v[0]  *dx + C.Grad.v[1]  *dy;
    double dp = C.Grad.p[0]  *dx + C.Grad.p[1]  *dy;

    W.rho += C.limiter_rho * dr;
    W.u   += C.limiter_u   * du;
    W.v   += C.limiter_v   * dv;
    W.p   += C.limiter_p   * dp;

    return W;
}




void compute_Lax_Friedrichs_flux_MUSCL( face2D &face,cellist_2D &C, flux &result, material &mat)
{
    int L = face.owner;
    int R = face.neighbour;

    cell2D &CL = C.cell_list[L];
    cell2D &CR = C.cell_list[R];

    const point2D &xf = face.face_centroid;


    primitive WL = reconstruct_face_primitive_limited(CL, xf);
    primitive WR = reconstruct_face_primitive_limited(CR, xf);
 

    conserved UL, UR;
    compute_conserved_cell_center(WL, UL, mat);
    compute_conserved_cell_center(WR, UR, mat);

    flux FL, FR;
    compute_normal_flux(face, FL, WL, mat);
    compute_normal_flux(face, FR, WR, mat);

    double nx = face.n.pos[0];
    double ny = face.n.pos[1];

    double unL = WL.u * nx + WL.v * ny;
    double unR = WR.u * nx + WR.v * ny;

    double lambda = compute_lambda(unL, unR, WL.p, WR.p,WL.rho, WR.rho,mat.gamma);

    for (int p = 0; p < 4; ++p)
    {
        result.F[p] = 0.5 * (FL.F[p] + FR.F[p]) - 0.5 * lambda * (UR.U[p] - UL.U[p]);
    }
    
}



void reset_residuals(cellist_2D &C)
{
    
    for (size_t c = 0; c < C.cell_list.size(); ++c)
    {
        for (int p = 0; p < 4; ++p)
            C.cell_list[c].residual[p] = 0.0;
    }
}


void compute_residual(cellist_2D &C, boundary_marker_list &boundary, facelist_2D &F, material &mat, simparam &sim)
{

    reset_residuals(C);
    // loop over faces

    for (size_t c = 0; c < C.cell_list.size(); ++c)
    {
        C.cell_list[c].Grad = {};   // zero explicitly
        compute_Gauss_Green_gradient(c, F, C);
    }

    for (size_t c = 0; c < C.cell_list.size(); ++c)
    {
        compute_venkat_limiter_cell(c, F, C);
    }

    for (size_t i = 0; i < F.face_list.size(); i++)
    {

        face2D &current_face = F.face_list[i];

        flux temp_result;

        int owner_index = current_face.owner;

        int neighbour_index = current_face.neighbour;

       

        if (neighbour_index != -1) // non boundary cells
        {
            grad temp;

            compute_Lax_Friedrichs_flux_MUSCL(current_face, C, temp_result, mat);

             
            for (int p = 0; p < 4; p++)
            {
                 C.cell_list[owner_index].residual[p] += temp_result.F[p] * current_face.len;
                 C.cell_list[neighbour_index].residual[p] -= temp_result.F[p] * current_face.len;
            }
        }
              
    }

    // compute the time step for all cells
    for (size_t i = 0; i < C.cell_list.size(); ++i)
    {
        C.cell_list[i].dt = compute_cell_time_step(C.cell_list[i], F, sim.CFL, mat);
    }
    
}

void compute_conserved_cell_center(primitive &var, conserved &cons, material &mat)
{

    cons.U[0] = var.rho;
    cons.U[1] = var.u * var.rho;
    cons.U[2] = var.v * var.rho;

    double vmag_squared = var.u * var.u + var.v * var.v;

    cons.U[3] = (var.p / (mat.gamma - 1)) + 0.5 * var.rho * vmag_squared;
}

void compute_fluxes_cell_center(flux_2D &flux, primitive &var, material &mat)
{

    double vmag_squared = var.u * var.u + var.v * var.v;

    flux.F[0] = var.rho * var.u;
    flux.F[1] = var.rho * var.u * var.u + var.p;
    flux.F[2] = var.rho * var.u * var.v;
    flux.F[3] = var.u * (var.p + (var.p / (mat.gamma - 1) + 0.5 * var.rho * vmag_squared));

    flux.G[0] = var.rho * var.v;
    flux.G[1] = var.rho * var.u * var.v;
    flux.G[2] = var.rho * var.v * var.v + var.p;
    flux.G[3] = var.v * (var.p + (var.p / (mat.gamma - 1) + 0.5 * var.rho * vmag_squared));
}


void compute_normal_flux(face2D &current_face, flux &Fn, primitive &var, material &mat)
{

    double nx = current_face.n.pos[0];
    double ny = current_face.n.pos[1];

    double un = var.u * nx + var.v * ny;

    double vmag_squared = var.u * var.u + var.v * var.v;

    Fn.F[0] = var.rho * un;
    Fn.F[1] = var.rho * var.u * un + var.p * nx;
    Fn.F[2] = var.rho * var.v * un + var.p * ny;
    Fn.F[3] = un * (var.p + (var.p / (mat.gamma - 1) + 0.5 * var.rho * vmag_squared));

}


void compute_Lax_Friedrichs_flux(face2D &current_face, cellist_2D &C, flux &result, material &mat)
{

    int Lcell_index, Rcell_index;

    Lcell_index = current_face.owner;
    Rcell_index = current_face.neighbour;

    cell2D Lcell = C.cell_list[Lcell_index];
    cell2D Rcell = C.cell_list[Rcell_index];

    compute_conserved_cell_center(Lcell.prim, Lcell.cons, mat);
    compute_conserved_cell_center(Rcell.prim, Rcell.cons, mat);

    flux Fn_L;
    flux Fn_R;

    compute_normal_flux(current_face, Fn_L, Lcell.prim, mat);
    compute_normal_flux(current_face, Fn_R, Rcell.prim, mat);

    double temp_lambda = 0;

    double nx = current_face.n.pos[0];
    double ny = current_face.n.pos[1];

    double un_L = Lcell.prim.u * nx + Lcell.prim.v * ny;

    double un_R = Rcell.prim.u * nx + Rcell.prim.v * ny;

    temp_lambda = compute_lambda(un_L, un_R, Lcell.prim.p, Rcell.prim.p, Lcell.prim.rho, Rcell.prim.rho, mat.gamma);

    for (int p = 0; p < 4; p++)
    {
        result.F[p] = 0.5 * (Fn_L.F[p] + Fn_R.F[p]) - 0.5 * temp_lambda * (Rcell.cons.U[p] - Lcell.cons.U[p]);
    }
     
}



void compute_Lax_Friedrichs_flux_freestream(face2D &current_face, cellist_2D &C, cell2D &ghost_neighbour, flux &result, material &mat, freestream &free_stream)
{

    int Lcell_index;

    Lcell_index = current_face.owner;

    cell2D Lcell = C.cell_list[Lcell_index];

    cell2D Rcell;

    double nx = current_face.n.pos[0];
    double ny = current_face.n.pos[1];

    double un_L = Lcell.prim.u * nx + Lcell.prim.v * ny;


    ghost_neighbour.prim.rho = free_stream.rho_inf;
    ghost_neighbour.prim.u =  free_stream.u_inf;
    ghost_neighbour.prim.v =  free_stream.v_inf;
    ghost_neighbour.prim.p =  free_stream.p_inf;


    compute_conserved_cell_center(Lcell.prim, Lcell.cons, mat);

    compute_conserved_cell_center(ghost_neighbour.prim, ghost_neighbour.cons, mat);

    flux Fn_L;
    flux Fn_R;


    compute_normal_flux(current_face, Fn_L, Lcell.prim, mat);
    compute_normal_flux(current_face, Fn_R, ghost_neighbour.prim, mat);

    Rcell = ghost_neighbour;

    double temp_lambda = 0;

    double un_R = Rcell.prim.u * nx + Rcell.prim.v * ny;

    temp_lambda = compute_lambda(un_L, un_R, Lcell.prim.p, Rcell.prim.p, Lcell.prim.rho, Rcell.prim.rho, mat.gamma);

    for (int p = 0; p < 4; p++)
    {

        result.F[p] = 0.5 * (Fn_L.F[p] + Fn_R.F[p]) - 0.5 * temp_lambda * (Rcell.cons.U[p] - Lcell.cons.U[p]);
    }

}



void compute_Lax_Friedrichs_flux_Euler_wall(face2D &current_face, cellist_2D &C, cell2D &ghost_neighbour, flux &result, material &mat, freestream &free_stream)
{

    int Lcell_index;

    Lcell_index = current_face.owner;

    cell2D Lcell = C.cell_list[Lcell_index];

    cell2D Rcell;

    compute_conserved_cell_center(Lcell.prim, Lcell.cons, mat);

    double rho, u, v, p, E;

    double nx = current_face.n.pos[0];
    double ny = current_face.n.pos[1];

    double un = Lcell.prim.u * nx + Lcell.prim.v * ny;

    rho = Lcell.cons.U[0];
    u = Lcell.cons.U[1] / rho;
    v = Lcell.cons.U[2] / rho;
    E = Lcell.cons.U[3] / rho;
    p = Lcell.prim.p;

    double ug = u - 2.0 * un * nx;
    double vg = v - 2.0 * un * ny;

    ghost_neighbour.cons.U[0] = rho;
    ghost_neighbour.cons.U[1] = rho * ug;
    ghost_neighbour.cons.U[2] = rho * vg;
    ghost_neighbour.cons.U[3] = p / (mat.gamma - 1.0) + 0.5 * rho * (ug * ug + vg * vg);

    ghost_neighbour.prim.rho = rho;
    ghost_neighbour.prim.u = ghost_neighbour.cons.U[1] / rho;
    ghost_neighbour.prim.v = ghost_neighbour.cons.U[2] / rho;
    ghost_neighbour.prim.p = (mat.gamma - 1.0) * (ghost_neighbour.cons.U[3] - 0.5 * rho * (ghost_neighbour.prim.u * ghost_neighbour.prim.u + ghost_neighbour.prim.v * ghost_neighbour.prim.v));

    flux Fn_L;
    flux Fn_R;


    compute_normal_flux(current_face, Fn_L, Lcell.prim, mat);
    compute_normal_flux(current_face, Fn_R, ghost_neighbour.prim, mat);

    Rcell = ghost_neighbour;

    double temp_lambda = 0;

    double un_L = Lcell.prim.u * nx + Lcell.prim.v * ny;

    double un_R = Rcell.prim.u * nx + Rcell.prim.v * ny;

    temp_lambda = compute_lambda(un_L, un_R, Lcell.prim.p, Rcell.prim.p, Lcell.prim.rho, Rcell.prim.rho, mat.gamma);

    for (int p = 0; p < 4; p++)
    {
        result.F[p] = 0.5 * (Fn_L.F[p] + Fn_R.F[p]) - 0.5 * temp_lambda * (Rcell.cons.U[p] - Lcell.cons.U[p]);
    }

}



double compute_cell_time_step(cell2D &C, facelist_2D &F, double CFL, material &mat)
{
    double cell_dt = 0;

    double denom = 0;

    for (int f = 0; f < 3; f++)
    {

        int face_id = C.face_id_indices[f];
        double un = C.prim.u * F.face_list[face_id].n.pos[0] + C.prim.v * F.face_list[face_id].n.pos[1];
        double a = sqrt(mat.gamma * C.prim.p / C.prim.rho);
        denom += (std::abs(un) + a) * F.face_list[face_id].len;

    }

    cell_dt = CFL * C.cell_area / denom;


    return cell_dt;
}

void solution_update(cellist_2D &C)
{
    for (size_t i = 0; i < C.cell_list.size(); ++i)
    {
        for (int p = 0; p < 4; p++)
        {
            C.cell_list[i].cons.U[p] = C.cell_list[i].cons.U[p] + (-C.cell_list[i].residual[p] / C.cell_list[i].cell_area) * C.cell_list[i].dt;
        }
    }
}



void apply_freestream(boundary_marker_list &boundary, cellist_2D &C, facelist_2D &F, freestream &free_stream, material &mat, int marker_index)
{

  for (size_t f = 0; f < boundary.marker_list[marker_index].facelist.size(); f++)
    {

        int fid = boundary.marker_list[marker_index].facelist[f];

        face2D &current_face = F.face_list[fid];

        int owner_index = current_face.owner;

        int neighbour_index = current_face.neighbour;

        cell2D ghost_neighbour;

        flux temp_result;

        compute_Lax_Friedrichs_flux_freestream(current_face, C, ghost_neighbour, temp_result, mat, free_stream);

        for (int p = 0; p < 4; p++)
        {
            C.cell_list[owner_index].residual[p] += temp_result.F[p] * current_face.len; // IMPORTANT, update the residual only for the owner cell
        }
    }

}


void apply_Euler_wall(boundary_marker_list &boundary, cellist_2D &C, facelist_2D &F, freestream &free_stream, material &mat,int marker_index)
{
   for (size_t f = 0; f < boundary.marker_list[marker_index].facelist.size(); f++)
    {

        int fid = boundary.marker_list[marker_index].facelist[f];

        face2D &current_face = F.face_list[fid];

        int owner_index = current_face.owner;

        int neighbour_index = current_face.neighbour;

        cell2D ghost_neighbour;

        flux temp_result;

        compute_Lax_Friedrichs_flux_Euler_wall(current_face, C, ghost_neighbour, temp_result, mat, free_stream);

        for (int p = 0; p < 4; p++)
        {
            C.cell_list[owner_index].residual[p] += temp_result.F[p] * current_face.len; // IMPORTANT, update the residual only for the owner cell
        }
    }
    
}


void apply_boundary_conditions(boundary_marker_list &boundary, cellist_2D &C, facelist_2D &F, freestream &free_stream, material &mat)
{

    for(int m=0;m<boundary.nmark;m++)
    {
        if(boundary.marker_list[m].btype==0)
        {
            apply_freestream(boundary,C,F,free_stream,mat,m);
        }
        else
        {
            apply_Euler_wall(boundary,C,F,free_stream,mat,m);
        }
    }

}

double compute_lambda(double un_L, double un_R, double P_L, double P_R, double rho_L, double rho_R, double gamma)
{

    double a_L = sqrt(gamma * P_L / rho_L);
    double a_R = sqrt(gamma * P_R / rho_R);

    return std::max(fabs(un_L) + a_L, fabs(un_R) + a_R);
}


std::array<double, 4> compute_L2_residual(cellist_2D &C)
{
    std::array<double, 4> R = {0, 0, 0, 0};


    for (auto &cell : C.cell_list)
    {
        double invA = 1.0 / cell.cell_area;

        for (int k = 0; k < 4; k++)
            R[k] += std::pow(cell.residual[k] * invA, 2);
    }

    for (int k = 0; k < 4; k++)
        R[k] = std::sqrt(R[k] / C.cell_list.size());

    return R;
}





