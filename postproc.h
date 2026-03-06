void write_vtk_2d(const std::string& filename, const std::vector<cell2D>& cells);
inline bool same_point(const point2D& a, const point2D& b);
void compute_pressure_coefficient(std::string filename ,boundary_marker_list &boundary,facelist_2D &F ,
    cellist_2D &clist,material &mat,freestream &free_stream, int marker_index);



inline bool same_point(const point2D& a, const point2D& b)
{
    const double eps = 1e-12;
    return std::abs(a.pos[0] - b.pos[0]) < eps &&
           std::abs(a.pos[1] - b.pos[1]) < eps;
}



void write_vtk_2d(const std::string& filename, const std::vector<cell2D>& cells)
{
    std::ofstream file(filename);
    if (!file) {
        throw std::runtime_error("Cannot open VTK file");
    }

    std::vector<point2D> points;                 // global unique points
    std::vector<std::array<int,3>> connectivity; // triangle connectivity

    // ---------- Build point list ----------
    for (const auto& c : cells) {
        std::array<int,3> conn;

        for (int k = 0; k < 3; ++k) {
            const point2D& p = c.tri.v[k];

            int id = -1;
            for (int i = 0; i < (int)points.size(); ++i) {
                if (same_point(points[i], p)) {
                    id = i;
                    break;
                }
            }

            if (id == -1) {
                id = points.size();
                points.push_back(p);
            }

            conn[k] = id;
        }

        connectivity.push_back(conn);
    }

    const int nPoints = points.size();
    const int nCells  = cells.size();

    // ---------- VTK header ----------
    file << "# vtk DataFile Version 3.0\n";
    file << "2D Euler FV solution\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";

    // ---------- Points ----------
    file << "POINTS " << nPoints << " float\n";
    for (const auto& p : points) {
        file << p.pos[0] << " "
             << p.pos[1] << " "
             << "0.0\n";
    }

    // ---------- Cells ----------
    file << "CELLS " << nCells << " " << 4 * nCells << "\n";
    for (const auto& conn : connectivity) {
        file << "3 "
             << conn[0] << " "
             << conn[1] << " "
             << conn[2] << "\n";
    }

    // ---------- Cell types ----------
    file << "CELL_TYPES " << nCells << "\n";
    for (int i = 0; i < nCells; ++i) {
        file << "5\n"; // VTK_TRIANGLE
    }

    // ---------- Cell data ----------
    file << "CELL_DATA " << nCells << "\n";

    // Density
    file << "SCALARS density float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (const auto& c : cells) {
        file << c.prim.rho << "\n";
    }

    // Velocity
    file << "VECTORS velocity float\n";
    for (const auto& c : cells) {
        file << c.prim.u << " "
             << c.prim.v << " 0.0\n";
    }

    // Pressure
    file << "SCALARS pressure float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (const auto& c : cells) {
        file << c.prim.p << "\n";
    }

    // Mach number
file << "SCALARS Mach float 1\n";
file << "LOOKUP_TABLE default\n";

const double gamma = 1.4;

for (const auto& c : cells) {
    double u = c.prim.u;
    double v = c.prim.v;
    double rho = c.prim.rho;
    double p = c.prim.p;

    double a = std::sqrt(gamma * p / rho);
    double V = std::sqrt(u*u + v*v);

    double M = V / a;

    file << M << "\n";
}


    file.close();
}



void compute_pressure_coefficient(std::string filename ,boundary_marker_list &boundary,facelist_2D &F ,
    cellist_2D &clist,material &mat,freestream &free_stream, int marker_index)
{

    std::fstream cp;

    double Fx=0;
    double Fy=0;

    double a_inf = std::sqrt(mat.gamma * mat.R * free_stream.T_inf);
    double V_inf = free_stream.Mach * a_inf;
    double rho  = free_stream.rho_inf;

    cp.open(filename,std::ios::out);

    double Cl=0;
    double Cd=0;

    for (size_t f = 0; f < boundary.marker_list[marker_index].facelist.size(); f++)
    {

    
        int fid = boundary.marker_list[marker_index].facelist[f];

        face2D current_face = F.face_list[fid];

        int owner_index = current_face.owner;

        double temp_cp= (clist.cell_list[owner_index].prim.p-free_stream.p_inf)/(0.5*free_stream.rho_inf*V_inf*V_inf);

        cp<<current_face.face_centroid.pos[0]<<","<<temp_cp<<std::endl;

        double nx = current_face.n.pos[0];
        double ny = current_face.n.pos[1];
        

        Fx+=  -(clist.cell_list[owner_index].prim.p)*current_face.len*nx;
        Fy+=  -(clist.cell_list[owner_index].prim.p)*current_face.len*ny;

        Cl +=  temp_cp * ny * current_face.len;
        Cd +=  temp_cp * nx * current_face.len;


    }
    double ca = cos(free_stream.AoA);
    double sa = sin(free_stream.AoA);

    double denom = 0.5*rho*V_inf*V_inf  ; 

    boundary.marker_list[marker_index].coeffs[0] = Cd;
    boundary.marker_list[marker_index].coeffs[1] = Cl;

    //std::cout << "Cl= " << Cl<<" , "
        //   << "Cd= "
        //   <<Cd
        //   << std::endl;




   cp.close();

}