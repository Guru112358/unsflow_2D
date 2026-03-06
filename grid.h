#include <cstdint>

struct grad
{
    std::array<double, 2>  rho; 
    std::array<double, 2>  u;
    std::array<double, 2>  v;
    std::array<double, 2>  p;   
};



struct cell2D
{
    int c_id;
    grad Grad;
    triangle tri;
    point2D centroid;
    std::array<int,3>  point_id;
    double cell_area; 
    int face_id_indices[3];  //facelist id numbers for 3 faces
    std::array<int,3> neighbour_cell_number; //element ids of neighbouring cells
    conserved cons;
    primitive prim;
    std::array<double, 4> residual;
    double dt;  //time step for each cell
    double limiter_rho = 1.0;
    double limiter_u   = 1.0;
    double limiter_v   = 1.0;
    double limiter_p   = 1.0;
};


struct cellist_2D
{
    std::vector<cell2D> cell_list;
};


struct face2D
{
    int f_id;
    std::array<int,2> vertex;  //contains the two vertex ids
    point2D face_centroid;
    point2D n;  //normal vector
    double len;  //face length
    int owner=-1;
    int neighbour=-1;
    flux_2D flux;
};

struct facelist_2D
{
    std::vector<face2D> face_list;
};


struct boundary_marker
{
    using intarray2 = std::array<int, 2>;
    std::string marker_tag;
    size_t n_elems;
    std::vector<intarray2> marker_nodes;
    std::vector<int> facelist;
    int btype ; // (0 for far field 1 for Euler wall)
    std::array<double,2> coeffs ; //cd - 0 index,  cl  - 1 index
    
};

struct boundary_marker_list
{
    int nmark;
    std::vector<boundary_marker> marker_list;
};

