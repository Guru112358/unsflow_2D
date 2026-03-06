
struct material
{
    double R;
    double gamma;
};


struct freestream
{
    double E_inf;
    double p_inf;
    double u_inf;
    double v_inf;
    double rho_inf;
    double T_inf;
    double Mach;
    double AoA;

};

freestream compute_freestream(freestream &fs,material &mat)

{

    fs.rho_inf = fs.p_inf / (mat.R * fs.T_inf);

    // Speed of sound
    double a_inf = std::sqrt(mat.gamma * mat.R * fs.T_inf);

    // Velocity magnitude
    double V_inf = fs.Mach * a_inf;

    // Angle of attack (radians)
    double alpha = fs.AoA * M_PI / 180.0;

    // Velocity components
    fs.u_inf =   V_inf * std::cos(alpha);
    fs.v_inf =   V_inf * std::sin(alpha);

    fs.E_inf = fs.p_inf / ((mat.gamma - 1.0) * fs.rho_inf) + 0.5 * (fs.u_inf * fs.u_inf + fs.v_inf * fs.v_inf);

    return fs;
    
}
