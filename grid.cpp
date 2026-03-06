#include "grid.h"
#include <unordered_map>

void read_grid(cellist_2D &clist, boundary_marker_list &boundary, std::string filename);
void compute_centroid_and_grid_volume_2D(cellist_2D &C);
void form_faces(cellist_2D &C, facelist_2D &F);
void find_neighbours(cellist_2D &C, facelist_2D &F);
void set_boundary_faces(boundary_marker_list &boundary, facelist_2D &F);
inline uint64_t edge_key(int a, int b);






inline uint64_t edge_key(int a, int b)
{
    if (a > b) std::swap(a, b);
    return (uint64_t(a) << 32) | uint64_t(b);
}




void read_grid(cellist_2D &clist, boundary_marker_list &boundary, std::string filename)
{

    int NELEM = 0;
    int NPOIN = 0;
    int NMARK = 0;
    std::ifstream in(filename);
    std::string line;

    while (std::getline(in, line))
    {
        if (line.find("NELEM=") != std::string::npos)
        {
            std::stringstream ss(line.substr(line.find("=") + 1));
            ss >> NELEM;
            break;
        }
    }

    clist.cell_list.resize(NELEM);

    for (int i = 0; i < NELEM; i++)
    {
        std::getline(in, line);
        std::stringstream ss(line);
        int type, id;
        int v0, v1, v2;
        ss >> type >> v0 >> v1 >> v2 >> id;

        clist.cell_list[i].point_id[0] = v0;
        clist.cell_list[i].point_id[1] = v1;
        clist.cell_list[i].point_id[2] = v2;

        clist.cell_list[i].c_id = id;

    }

    while (std::getline(in, line))
    {
        if (line.find("NPOIN") != std::string::npos)
        {
            std::stringstream ss(line.substr(line.find("=") + 1));
            ss >> NPOIN;
            break;
        }
    }

    std::vector<point2D> ps;
    ps.resize(NPOIN);

    for (int i = 0; i < NPOIN; i++)
    {
        std::getline(in, line);
        std::stringstream ss(line);
        int p;
        double temp_x, temp_y;

        ss >> temp_x >> temp_y >> p;
        ps[i].pos[0] = temp_x;
        ps[i].pos[1] = temp_y;
    }

    // reconstitute grid

    for (int i = 0; i < NELEM; i++)
    {

        clist.cell_list[i].tri.v[0] = {ps[clist.cell_list[i].point_id[0]].pos[0], ps[clist.cell_list[i].point_id[0]].pos[1]};

        clist.cell_list[i].tri.v[1] = {ps[clist.cell_list[i].point_id[1]].pos[0], ps[clist.cell_list[i].point_id[1]].pos[1]};

        clist.cell_list[i].tri.v[2] = {ps[clist.cell_list[i].point_id[2]].pos[0], ps[clist.cell_list[i].point_id[2]].pos[1]};
    }

    while (std::getline(in, line))
    {
        if (line.find("NMARK=") != std::string::npos)
        {
            std::stringstream ss(line.substr(line.find("=") + 1));
            ss >> NMARK;
            break;
        }
    }
    boundary.marker_list.resize(NMARK);
    boundary.nmark = NMARK;

    // Loop through all NMARK boundary blocks
    for (int m = 0; m < NMARK; ++m)
    {
        std::string line;
        std::string tag_name;
        int num_elements;

        std::getline(in, line);
        size_t value_start = line.find("=") + 1;
        tag_name = line.substr(value_start);
        tag_name.erase(0, tag_name.find_first_not_of(" \t\n\r"));

        std::getline(in, line);
        value_start = line.find("=") + 1;
        std::stringstream ss_elems(line.substr(value_start));
        ss_elems >> num_elements;

        boundary.marker_list[m].marker_tag = tag_name;
        boundary.marker_list[m].marker_nodes.resize(num_elements);

        for (int i = 0; i < num_elements; ++i)
        {
            std::getline(in, line);
            std::stringstream ss_edge(line);

            int element_type, node_a, node_b;
            
            ss_edge >> element_type >> node_a >> node_b;

            if (node_a > node_b)
                std::swap(node_a, node_b); 


            boundary.marker_list[m].marker_nodes[i] = {node_a, node_b};
        }
    }
}

void compute_centroid_and_grid_volume_2D(cellist_2D &C)
{

    // computing centroid of the triangle
    double one_third = 1.0 / 3;
    for (size_t i = 0; i < C.cell_list.size(); ++i)
    {
        C.cell_list[i].centroid.pos[0] = one_third * (C.cell_list[i].tri.v[0].pos[0] + C.cell_list[i].tri.v[1].pos[0] + C.cell_list[i].tri.v[2].pos[0]);
        C.cell_list[i].centroid.pos[1] = one_third * (C.cell_list[i].tri.v[0].pos[1] + C.cell_list[i].tri.v[1].pos[1] + C.cell_list[i].tri.v[2].pos[1]);
    }

    // computing the area using the expression for the area of a triangle

    for (size_t i = 0; i < C.cell_list.size(); ++i)
    {
        C.cell_list[i].cell_area =
            0.5 *
            (C.cell_list[i].tri.v[0].pos[0] * (C.cell_list[i].tri.v[1].pos[1] - C.cell_list[i].tri.v[2].pos[1]) +
             C.cell_list[i].tri.v[1].pos[0] * (C.cell_list[i].tri.v[2].pos[1] - C.cell_list[i].tri.v[0].pos[1]) +
             C.cell_list[i].tri.v[2].pos[0] * (C.cell_list[i].tri.v[0].pos[1] - C.cell_list[i].tri.v[1].pos[1])

            );

    }
}



void form_faces(cellist_2D &C, facelist_2D &F)
{
    F.face_list.clear();

    size_t Ncells = C.cell_list.size();

    std::unordered_map<uint64_t, int> edge2face;

    edge2face.reserve(3 * Ncells);

    for (size_t i = 0; i < Ncells; ++i)
    {
        auto &cell = C.cell_list[i];

        int vids[3] = {
            cell.point_id[0],
            cell.point_id[1],
            cell.point_id[2]
        };

        for (int f = 0; f < 3; ++f)
        {
            int v0 = vids[f];
            int v1 = vids[(f + 1) % 3];

            uint64_t key = edge_key(v0, v1);

            auto it = edge2face.find(key);

            if (it == edge2face.end())
            {

                face2D fc;
                fc.vertex[0] = std::min(v0, v1);
                fc.vertex[1] = std::max(v0, v1);
                fc.owner     = i;
                fc.neighbour = -1;

                const point2D &p0 = cell.tri.v[f];
                const point2D &p1 = cell.tri.v[(f + 1) % 3];

                fc.face_centroid.pos[0] = 0.5 * (p0.pos[0] + p1.pos[0]);
                fc.face_centroid.pos[1] = 0.5 * (p0.pos[1] + p1.pos[1]);

                double dx = p1.pos[0] - p0.pos[0];
                double dy = p1.pos[1] - p0.pos[1];
                fc.len = std::sqrt(dx * dx + dy * dy);

                fc.n.pos[0] =  dy / fc.len;
                fc.n.pos[1] = -dx / fc.len;

                double cx = fc.face_centroid.pos[0] - cell.centroid.pos[0];
                double cy = fc.face_centroid.pos[1] - cell.centroid.pos[1];

                if (fc.n.pos[0] * cx + fc.n.pos[1] * cy < 0.0)
                {
                    fc.n.pos[0] *= -1.0;
                    fc.n.pos[1] *= -1.0;
                }

                int face_id = (int)F.face_list.size();
                F.face_list.push_back(fc);

                edge2face[key] = face_id;
                cell.face_id_indices[f] = face_id;
            }
            else
            {
                int face_id = it->second;
                face2D &fc = F.face_list[face_id];

                fc.neighbour = i;
                cell.face_id_indices[f] = face_id;
            }
        }
    }
    
}




void find_neighbours(cellist_2D &C, facelist_2D &F)
{
     int Ncells = C.cell_list.size();

     //set all  neigbours to -1 initially
    for (int i = 0; i < Ncells; ++i)
    {
        for (int f = 0; f < 3; ++f)
        {
            C.cell_list[i].neighbour_cell_number[f] = -1;
        }
    }
    
    for (size_t face_id = 0; face_id < F.face_list.size(); ++face_id)
    {
        const face2D &fc = F.face_list[face_id];

        if (fc.owner >= 0 && fc.owner < Ncells)
        {
            for (int lf = 0; lf < 3; ++lf)
            {
                if (C.cell_list[fc.owner].face_id_indices[lf] == (int)face_id)
                {
                    C.cell_list[fc.owner].neighbour_cell_number[lf] = fc.neighbour;
                    break;
                }
            }
        }

        if (fc.neighbour >= 0 && fc.neighbour < Ncells)
        {
            for (int lf = 0; lf < 3; ++lf)
            {
                if (C.cell_list[fc.neighbour].face_id_indices[lf] == (int)face_id)
                {
                    C.cell_list[fc.neighbour].neighbour_cell_number[lf] = fc.owner;
                    break;
                }
            }
        }
    }
}


// this has bad time complexity but in pracice number of boundary faces are not that big 
//that I care to change this

void set_boundary_faces(boundary_marker_list &boundary, facelist_2D &F)
{

    for (int m = 0; m < boundary.nmark; m++)
    {

        for (size_t i = 0; i < boundary.marker_list[m].marker_nodes.size(); i++)
        {

            int a = boundary.marker_list[m].marker_nodes[i][0];
            int b = boundary.marker_list[m].marker_nodes[i][1];

            // now loop over faces to find the correct face id
            for (size_t f = 0; f < F.face_list.size(); f++)
            {

                face2D current_face = F.face_list[f];

                if (a == current_face.vertex[0] && b == current_face.vertex[1])
                {
                    boundary.marker_list[m].facelist.push_back(f);

                }
            }
        }
    }
}
