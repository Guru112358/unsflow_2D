

struct conserved
{
std::array<double, 4> U;
};


struct flux_2D
{

 std::array<double, 4> F;
 std::array<double, 4> G;  

};


struct flux
{

 std::array<double, 4> F;  

};

struct primitive
{
    double rho;
    double u;
    double v;
    double p;
};