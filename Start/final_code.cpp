#include <algorithm>
#include <chrono>
using namespace std::chrono;
#include <limits>
#include <iostream>
#include <tuple>
#include <list>

#include <random>
static std::default_random_engine engine(10); // random seed = 10
static std::uniform_real_distribution<double> uniform(0, 1);

#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
 
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
 
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define PI 3.14159

double sqr(double x) {return x * x;};

void boxMuller(double stdev, double &x, double &y) {
    double r1 = uniform(engine);
    double r2 = uniform(engine);
    x = sqrt(-2 * log(r1)) * cos(2 * PI * r2) * stdev;
    y = sqrt(-2 * log(r1)) * sin(2 * PI * r2) * stdev;

}
 
class Vector {
public:
    explicit Vector(double x = 0, double y = 0, double z = 0) {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }
    double norm2() const {
        return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
    }
    double norm() const {
        return sqrt(norm2());
    }
    void normalize() {
        double n = norm();
        data[0] /= n;
        data[1] /= n;
        data[2] /= n;
    }
    double operator[](int i) const { return data[i]; };
    double& operator[](int i) { return data[i]; };
    double data[3];
};
 
Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator-(const Vector& a) {
    return Vector(-a[0], -a[1], -a[2]);
}
Vector operator*(const double a, const Vector& b) {
    return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const double b) {
    return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator/(const Vector& a, const double b) {
    return Vector(a[0] / b, a[1] / b, a[2] / b);
}
Vector operator*(const Vector& a, const Vector& b) {
    return Vector(a[0]*b[0], a[1]*b[1], a[2]*b[2]);
}
double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

Vector random_cos(const Vector &N) {
    double r1 = uniform(engine);
    double r2 = uniform(engine);
    double x = cos(2 * PI * r1) * sqrt(1 - r2);
    double y = sin(2 * PI * r1) * sqrt(1 - r2);
    double z = sqrt(r2);

    Vector T1;
    double min_value_index_N = 0;
    double min_value = std::numeric_limits<double>::max();
    for (int i = 0; i < 3; i++) {
        if (abs(N[i]) < abs(N[min_value_index_N])) {
            min_value_index_N = i;
            min_value = abs(N[i]);
        }
    }

    T1[min_value_index_N] = 0.;
    int index_swap_1 = abs(min_value_index_N - 1);
    int index_swap_2 = abs(min_value_index_N - 2);
    double temp = T1[index_swap_1];
    T1[index_swap_1] = T1[index_swap_2];
    T1[index_swap_2] = -1 * temp;
    
    Vector T2 = cross(N, T1);
    Vector V = x * T1 + y * T2 + z * N;
    return V;
}

class Ray {
public:
    Ray(const Vector& O, const Vector &u): origin(O), u(u) {};
    Vector origin;
    Vector u;
};

class Object {
    public:
        Object(const Vector& albedo = Vector(1., 1., 1.), bool isMirror = false, bool isTransparent = false): albedo(albedo), isMirror(isMirror), isTransparent(isTransparent) {}
        virtual bool intersect(const Ray& r, Vector &P, Vector&N, double &t) const = 0;

        Vector albedo;
        bool isMirror;
        bool isTransparent;
};

// _________________________________________________________________________________________________________________________________________________
// Parts of this code was taken from: https://pastebin.com/CAgp9r15

#include <string>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <vector>

class BoundingBox {
public:
    Vector B_min;
    Vector B_max;

    explicit BoundingBox(Vector min = Vector(), Vector max = Vector()) {
        B_min = min;
        B_max = max;
    }

    bool bounding_box_intersects(const Ray& r, double& t) const {
        double t_min = -std::numeric_limits<double>::infinity();
        double t_max = std::numeric_limits<double>::infinity();
    
        for (int i = 0; i < 3; ++i) {
            double t0 = (B_min[i] - r.origin[i]) / r.u[i];
            double t1 = (B_max[i] - r.origin[i]) / r.u[i];
            if (t0 > t1) std::swap(t0, t1);
            t_min = std::max(t_min, t0);
            t_max = std::min(t_max, t1);
    
            if (t_max < t_min || t_max < 0) return false;
        }
    
        t = t_min;
        return true;
    }
    
};
 
class TriangleIndices {
public:
    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
    };
    int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
    int uvi, uvj, uvk;  // indices within the uv coordinates array
    int ni, nj, nk;  // indices within the normals array
    int group;       // face group
};
 
 
class TriangleMesh : public Object {
    // Vector color;
    // double refractive_index;
    // bool reflects;
public:
    class BVH {
    public:
        BoundingBox box;
        BVH* left, * right;
        int first, last;
        TriangleMesh* mesh;
    
        // Define a scale and translate function
    
        // Leal Koksal helped me with this function as there was a bug I did not manage to solve on my own
        BoundingBox computeBbox(int starting_triangle, int ending_triangle) {
            Vector vertex = this->mesh->scaling_factor * this->mesh->vertices[this->mesh->indices[starting_triangle].vtxk] + this->mesh->translation;
            BoundingBox initial_bbox = BoundingBox(vertex, vertex);

            for (int i=starting_triangle; i<ending_triangle; i++) {
                auto triangle_vertices = {
                    this->mesh->vertices[this->mesh->indices[i].vtxi], 
                    this->mesh->vertices[this->mesh->indices[i].vtxj],
                    this->mesh->vertices[this->mesh->indices[i].vtxk]
                };
                for (auto const& v : triangle_vertices) {
                    Vector vertex = this->mesh->scaling_factor * v + this->mesh->translation;

                    for (int k = 0; k < 3; ++k) {
                        initial_bbox.B_min[k] = std::min(initial_bbox.B_min[k], vertex[k]);
                        initial_bbox.B_max[k] = std::max(initial_bbox.B_max[k], vertex[k]);
                    }
                }
            }
            return initial_bbox;
        }
    
        void build_bvh(BVH* current, int first, int last) {
            current->first = first;
            current->last = last;
            current->box = computeBbox(first, last);
            // current->left = NULL;
            // current->right = NULL;
    
            if (last-first < 5) return;
    
            Vector diagonal = current->box.B_max - current->box.B_min;
    
            // Finding longest axis to make it efficient
            int axis = 2;
            if ((diagonal[0] >= diagonal[1]) && (diagonal[0] >= diagonal[2])) {
                axis = 0;
            }
            else {
                if ((diagonal[1] >= diagonal[0]) && (diagonal[1] >= diagonal[2])) {
                    axis = 1;
                }
            }
    
            double middle = current->box.B_min[axis] + 0.5 * diagonal[axis];
            int pivot = first;
            for (int i=first; i < last; i++) {
                // Barycentre
                double bary = (mesh->vertices[mesh->indices[i].vtxi][axis] + mesh->vertices[mesh->indices[i].vtxj][axis] + mesh->vertices[mesh->indices[i].vtxk][axis]) / 3;
                if (bary < middle) {
                    std::swap(mesh->indices[i], mesh->indices[pivot]);
                    pivot++;
                }
            }
            if ((pivot <=first) || (pivot >= last - 1)) return; // if pivot is too much to the left
    
            current->left = new BVH();
            current->right = new BVH();
            build_bvh(current->left, first, pivot);
            build_bvh(current->right, pivot, last);
        }
    
    };

    BVH bvh;
    double scaling_factor;
    Vector translation;
    ~TriangleMesh() {}
    TriangleMesh(double scaling_factor, Vector translation, Vector albedo, double isMirror = false, bool isTransparent = false) : Object(albedo, isMirror, isTransparent), scaling_factor(scaling_factor), translation(translation) {}; //, root(new Node) {};


    void readOBJ(const char* obj) {
        char matfile[255];
        char grp[255];
 
        FILE* f;
        f = fopen(obj, "r");
        int curGroup = -1;
        while (!feof(f)) {
            char line[255];
            if (!fgets(line, 255, f)) break;
 
            std::string linetrim(line);
            linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
            strcpy(line, linetrim.c_str());
 
            if (line[0] == 'u' && line[1] == 's') {
                sscanf(line, "usemtl %[^\n]\n", grp);
                curGroup++;
            }
 
            if (line[0] == 'v' && line[1] == ' ') {
                Vector vec;
 
                Vector col;
                if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
                    col[0] = std::min(1., std::max(0., col[0]));
                    col[1] = std::min(1., std::max(0., col[1]));
                    col[2] = std::min(1., std::max(0., col[2]));
 
                    vertices.push_back(vec);
                    vertexcolors.push_back(col);
 
                } else {
                    sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                    vertices.push_back(vec);
                }

                Vector transformed = scaling_factor * vec + translation;

                // expand the bounding box
                for (int k = 0; k < 3; ++k) {
                    bbox.B_min[k] = std::min(bbox.B_min[k], transformed[k]);
                    bbox.B_max[k] = std::max(bbox.B_max[k], transformed[k]);
                }
            }
            if (line[0] == 'v' && line[1] == 'n') {
                Vector vec;
                sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                normals.push_back(vec);
            }
            if (line[0] == 'v' && line[1] == 't') {
                Vector vec;
                sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
                uvs.push_back(vec);
            }
            if (line[0] == 'f') {
                TriangleIndices t;
                int i0, i1, i2, i3;
                int j0, j1, j2, j3;
                int k0, k1, k2, k3;
                int nn;
                t.group = curGroup;
 
                char* consumedline = line + 1;
                int offset;
 
                nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
                if (nn == 9) {
                    if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                    if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                    if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                    if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                    if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                    if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                    if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                    if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                    if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                    indices.push_back(t);
                } else {
                    nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                    if (nn == 6) {
                        if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                        if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                        if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                        if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                        if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                        if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                        indices.push_back(t);
                    } else {
                        nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                        if (nn == 3) {
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            indices.push_back(t);
                        } else {
                            nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                            if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                            if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                            indices.push_back(t);
                        }
                    }
                }
 
                consumedline = consumedline + offset;
 
                while (true) {
                    if (consumedline[0] == '\n') break;
                    if (consumedline[0] == '\0') break;
                    nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
                    TriangleIndices t2;
                    t2.group = curGroup;
                    if (nn == 3) {
                        if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                        if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                        if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                        if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                        if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                        if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                        if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                        if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                        if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
                        indices.push_back(t2);
                        consumedline = consumedline + offset;
                        i2 = i3;
                        j2 = j3;
                        k2 = k3;
                    } else {
                        nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                        if (nn == 2) {
                            if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                            if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                            if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                            if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                            if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                            if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                            consumedline = consumedline + offset;
                            i2 = i3;
                            j2 = j3;
                            indices.push_back(t2);
                        } else {
                            nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                            if (nn == 2) {
                                if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                                if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                                if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;                             
                                consumedline = consumedline + offset;
                                i2 = i3;
                                k2 = k3;
                                indices.push_back(t2);
                            } else {
                                nn = sscanf(consumedline, "%u%n", &i3, &offset);
                                if (nn == 1) {
                                    if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                    if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                    if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                    consumedline = consumedline + offset;
                                    i2 = i3;
                                    indices.push_back(t2);
                                } else {
                                    consumedline = consumedline + 1;
                                }
                            }
                        }
                    }
                }
 
            }
        }
        fclose(f);

        bvh.mesh = this;
        bvh.build_bvh(&bvh, 0, indices.size());
    }

    // Bounding Box version
    /*
    bool intersect(const Ray& r, Vector &P, Vector&N, double &t) const override {
        // bool intersect(const Ray& r, Vector &P, Vector &N, double &t) {
        t = std::numeric_limits<double>::max();
        double t_box;
        if (!bbox.bounding_box_intersects(r, t_box)) {
            return false;
        }

        bool has_intersection = false;

        for (int i = 0; i < indices.size(); i++) {
            const Vector& A = scaling_factor * vertices[indices[i].vtxi] + translation; // put scaling and translation in a separate function
            const Vector& B = scaling_factor * vertices[indices[i].vtxj] + translation;
            const Vector& C = scaling_factor * vertices[indices[i].vtxk] + translation;

            const Vector e1 = B - A;
            const Vector e2 = C - A;
            //Vector localN =  scaling_factor * normals[indices[i].vtxi];
            //localN.normalize();
            Vector localN = cross (e1, e2);
            //localN.normalize();
            // double invUN = 1. / dot(r.u, N);
            Vector AOcrossU = cross(A-r.origin, r.u);

            double localt = dot(A-r.origin, localN) / dot(r.u, localN); // implement the formula from the slides
            // if (localt < 0) continue;
            // if (localt > t) continue;

            double beta = dot(e2, AOcrossU) / dot(r.u, localN); // implement formula from slides
            // check condition < 0 || > 1

            double gamma = -dot(e1, AOcrossU) / dot(r.u, localN); // implement the formula from the slides
            //check condition <0 || > 1

            double alpha = 1 - beta - gamma; // implement the formula from the slides
            // check condition <0 

            if (alpha >= 0.0f && alpha <= 1.0f && beta >= 0.0f && beta <= 1.0f && gamma >= 0.0f && gamma <= 1.0f && localt >= 0.0f && localt < t) {
                t = localt;
                P = r.origin + t * r.u;
                N = localN;
                N.normalize();
                has_intersection = true;
            }
        }
        return has_intersection;
    }
    */

    // BVH version
    bool intersect(const Ray& r, Vector& P, Vector& N, double& t) const {
        bool result = false;
        t = std::numeric_limits<double>::max();
        double t_box;
        if (!bvh.box.bounding_box_intersects(r, t_box)) return false;

        std::list<const BVH*> nodes_to_visit;
        nodes_to_visit.push_back(&bvh);
        while(!nodes_to_visit.empty()) {
            const BVH* cur = nodes_to_visit.back();
            nodes_to_visit.pop_back();
            if (cur->left) { // not a leaf
                // double t_box;
                if (cur->left->box.bounding_box_intersects(r, t_box)) {
                    if (t_box < t) {
                        nodes_to_visit.push_back(cur->left);
                    }
                }
                // double t_box;
                if (cur->right->box.bounding_box_intersects(r, t_box)) {
                    if (t_box < t) {
                        nodes_to_visit.push_back(cur->right);
                    }  
                }
            }
            else {
                // copy all the code
                bool has_intersection = false;
                for (int i = cur->first; i < cur->last; i++) {
                    const Vector& A = scaling_factor * vertices[indices[i].vtxi] + translation; // put scaling and translation in a separate function
                    const Vector& B = scaling_factor * vertices[indices[i].vtxj] + translation;
                    const Vector& C = scaling_factor * vertices[indices[i].vtxk] + translation;

                    const Vector e1 = B - A;
                    const Vector e2 = C - A;
                    Vector localN = cross (e1, e2);
                    Vector AOcrossU = cross(A-r.origin, r.u);

                    double localt = dot(A-r.origin, localN) / dot(r.u, localN); // implement the formula from the slides
                    double beta = dot(e2, AOcrossU) / dot(r.u, localN); // implement formula from slides
                    double gamma = -dot(e1, AOcrossU) / dot(r.u, localN); // implement the formula from the slides
                    double alpha = 1 - beta - gamma; // implement the formula from the slides

                    if (alpha >= 0.0f && alpha <= 1.0f && beta >= 0.0f && beta <= 1.0f && gamma >= 0.0f && gamma <= 1.0f && localt >= 0.0f && localt < t) {
                        t = localt;
                        P = r.origin + t * r.u;
                        N = localN;
                        N.normalize();
                        has_intersection = true;
                    }
                    // return has_intersection;
                }
                if (has_intersection) {
                    result = true;
                    // return true;
                }
            }
        }
        return result;
    }
 
    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;
    BoundingBox bbox;
    
};

// _________________________________________________________________________________________________________________________________________________
 
class Sphere : public Object {
public:
    Sphere(const Vector& C, double R, Vector albedo, bool isTransparent, bool isMirror, bool isHollow) 
        : Object(albedo, isMirror, isTransparent), C(C), R(R), is_hollow(isHollow) {}

    Sphere(const Vector& C, double R, Vector albedo, bool isTransparent, bool isMirror) 
        : Object(albedo, isMirror, isTransparent), C(C), R(R), is_hollow(false) {}
    // Sphere(const Vector& C, double R, Vector albedo, bool isTransparent, bool isMirror, bool isHollow) : C(C), R(R), Object(albedo, isMirror, isTransparent), is_hollow(isHollow) {}; // , bool is_H) : C(C), R(R), albedo(albedo), transparent(is_T), mirror(is_M), is_hollow(is_H) {};
    // Sphere(const Vector& C, double R, Vector albedo, bool is_T, bool is_M) : C(C), R(R), Object(albedo, isMirror, isTransparent), is_hollow(false) {}; // albedo(albedo), transparent(is_T), mirror(is_M), is_hollow(false) {};
    Vector C;
    double R;
    // Vector albedo;
    // bool transparent;
    // bool mirror;
    bool is_hollow;

    bool intersect(const Ray& r, Vector &P, Vector &N, double &t) const override { // here, returned by reference but can do something else
        double delta = sqr(dot(r.u, r.origin-C)) - ((r.origin-C).norm2() - sqr(R)); // discriminant formula from slides
        if (delta < 0) {
            return false; // we only wants points in front of us?
        }
        double x = dot(r.u, C-r.origin);
        double t1 = x - sqrt(delta);
        double t2 = x + sqrt(delta);
        if (t2 < 0) {
            return false; // there is not intersection
        }
        // otherwise either t1 or t2 is in front of me
        if (t1 >= 0) {
            t = t1;
        }
        else {
            t = t2;
        }
        // if (t < 1e-4) {  // Ignore self-intersections
        //     return false;
        // }
        P = r.origin + t * r.u;
        N = (P-C) / (P-C).norm();
        
        return true;
    }
};

class Scene {
    public:
        // Scene(std::vector<Sphere> obj) : objects(obj) {};
        Scene(std::vector<Object *> obj) : objects(obj) {};

        // std::vector<Sphere> objects;
        std::vector<Object *> objects;
        // add albedo
        // add light color?

        void add(Sphere *s) {
            objects.push_back(s);
        }

        bool intersect(const Ray& r, Vector& P, Vector& N, double& t, int &id) {
            bool result = false;
            double closest_t = std::numeric_limits<double>::max(); // maximum double so any sphere is closer initially
            
            for (int i = 0; i < objects.size(); i++) {

                Vector localN, localP; // values for current local object
                double localt;

                if (objects[i]->intersect(r, localP, localN, localt)) {
                    result = true;
                    if (localt < closest_t) {
                        N = localN;
                        P = localP;
                        closest_t = localt;
                        id = i;
                    }
                }
            }
            t = closest_t;
            return result;
        }

        /*
        void intersect(const Ray& r, Vector& P, Vector& N) {
            double t = std::numeric_limits<double>::max(); // wrong 
            for (int i = 0; i < objects.size(); i++) {

                Vector localN, localP; // values for current local object
                double localt;
                objects[i].intersect(r, localP, localN);
            }
        }
        */
        bool isShadow(Scene& scene, const Vector& light_position, const Vector& P, const Vector& N, double& t) {
            Vector lightDir = light_position - P;
            double lightDist = lightDir.norm();
            lightDir.normalize();
        
            Ray shadowRay = Ray(P + 1e-4 * N, lightDir);
            Vector shadowP, shadowN;
            double shadowt;
            int shadow_id;
        
            if (scene.intersect(shadowRay, shadowP, shadowN, shadowt, shadow_id)) {
                double intersectDist = (shadowP - P).norm();
                if (intersectDist < lightDist) {
                    return true;
                }
            }
            return false;
        }


        Vector getColor(const Ray &ray, int bounce_number, Scene& scene, const Vector& light_position, const double& I) {
            if (bounce_number <= 0) {return Vector(0, 0, 0);} // ? // stop the iteration //else, decrease bounce_number
            // put nearly everything in it

            Vector P, N;
            double t; // maximum double so any sphere is closer initially;
            int id;
            Vector color;

            if (intersect(ray, P, N, t, id)) {
                if (objects[id]->isMirror) {
                    // Ray reflected_ray = ...;
                    Vector reflectDir = ray.u - 2 * dot(ray.u, N) * N;
                    Ray reflectedRay(P + 0.00001 * N, reflectDir);
                    return getColor(reflectedRay, bounce_number-1, scene, light_position, I);
                }

                if (objects[id]->isTransparent) {
                    // Transparence
                    double n1 = 1;
                    double n2 = 1.5;
                    if (dot(ray.u, N) > 0) {
                        std::swap(n1, n2);
                        N = -N;
                    }
                    /*
                    if (objects[id]->is_hollow) {
                        std::swap(n1, n2);
                        N = -N;
                    }
                    */
                    Vector lightDir = light_position - P;
                    double lightDist = lightDir.norm();
                    Vector albedo = objects[id]->albedo;
                    Vector omega_i = ray.u; // incoming ray;
                    Vector tangentialComponent = (n1 / n2) * (omega_i - dot(omega_i, N) * N); // formula from slides
                    double D = 1. - std::pow(n1/n2, 2)*(1. - std::pow(dot(omega_i, N), 2)); // formula from slides // make sure stuff in 
                    Vector normalComponent;
                    if (D <= 0. ) {
                        // consider it to be a mirror
                        // normalComponent = -sqrt(-D) * N;
                        Vector reflectDir = ray.u - (2 * dot(ray.u, N) * N);
                        Ray reflectedRay(P + 1e-4 * N, reflectDir);
                        return getColor(reflectedRay, bounce_number-1, scene, light_position, I);
                    }
                    else {
                        normalComponent = -sqrt(D)*N;
                        Vector transmittedDirection = normalComponent + tangentialComponent; //
                        transmittedDirection.normalize();

                        // Fresnel's law
                        double k0 = pow(n1 - n2, 2) / pow(n1 + n2, 2);
                        double R = k0 + (1 - k0) * pow(1 - abs(dot(N, omega_i)), 5);
                        double T = 1- R;

                        // Solution (Leal Koksal gave me a hint to help code this)
                        Vector reflectDir = ray.u - (2 * dot(ray.u, N) * N);
                        Ray reflectedRay(P + 1e-4 * N, reflectDir);

                        Ray refracted_ray = Ray(P - 1e-4 * N, transmittedDirection);

                        return R * getColor(reflectedRay, bounce_number-1, scene, light_position, I) + T * getColor(refracted_ray, bounce_number - 1, scene, light_position, I);

                        /*
                        // Solution with random point which is very noisy
                        double u = uniform(engine); // random number between 0 and 1
                        if (u < R) {
                            // launch a reflection ray
                            Vector reflectDir = ray.u - (2 * dot(ray.u, N) * N);
                            Ray reflectedRay(P + 1e-4 * N, reflectDir);
                            return getColor(reflectedRay, bounce_number-1, scene, light_position, I);
                        }
                        else {
                            // launch a refraction ray
                            Ray refracted_ray = Ray(P - 1e-4 * N, transmittedDirection);
                            return getColor(refracted_ray, bounce_number - 1, scene, light_position, I);
                        }
                        */

                        // Without Fresnel's law
                        // Ray transmittedRay(P - 1e-4 * N, transmittedDirection); // go inside object (make sure you don't interest with same border)
                        // Vector transmittedColor = getColor(transmittedRay, bounce_number-1, scene, light_position, I); //
                        // return transmittedColor;


                    }                    
                }

                else {
                    // handle diffuse surfaces
                    Vector Lo(0., 0., 0.);
                    // add direct lighting;
                    double visibility = 1.; // computes the visibility term by launching a ray towards the light source
        
                    if (isShadow(scene, light_position, P, N, t)) {
                        visibility = 0.;
                    }
                
                    Vector lightDir = light_position - P;
                    double lightDist = lightDir.norm();
                    Vector albedo = objects[id]->albedo;
                    Vector omega_i = lightDir/lightDir.norm();
                    double dot_product = std::max(dot(N, omega_i), 0.); // make sure not negative
                    Lo = (I/(4 * PI * lightDir.norm2())) * (albedo / PI) * visibility * dot_product;

                    // add indirect lighting
                    Ray randomRay = Ray(P, random_cos(N));; // randomly sample ray using random_cos
                    Lo = Lo + albedo * getColor(randomRay, bounce_number-1, scene, light_position, I);
                    
                    /*
                        Vector lightDir = light_position - P;
                        double lightDist = lightDir.norm();
                        Vector albedo = objects[id].albedo;
                        Vector omega_i = lightDir/lightDir.norm();
                        double dot_product = std::max(dot(N, omega_i), 0.); // make sure not negative
                        color = (I/(4 * PI * lightDir.norm2())) * (albedo / PI) * dot_product;

                        if (isShadow(scene, light_position, P, N, t)) {
                            color = Vector(0, 0, 0);
                        }
                    */
                    return Lo;
                }
            }
            return color;
        }
};

 
int main() {
    // The code related to the timer was found on stack overflow at the following link:
    // https://stackoverflow.com/questions/30702759/implementing-a-timer-in-c
    // Start a timer
    auto a = high_resolution_clock::now();
    /*
        Vector del1 = Vector(2, -2, 1);
        std::cout << "del1 = " << del1[0] << ", " << del1[1] << ", " << del1[2] << std::endl;
        Vector del2 = random_cos(del1);
        std::cout << "del2 = " << del2[0] << ", " << del2[1] << ", " << del2[2] << std::endl;
    */
    int W = 256; //512;
    int H = 256; // 512;
    Vector camera_origin(0, 0, 55);
    double fov = 60 * PI / 180; // make sure radians not degrees

    Object* S_mirror = new Sphere(Vector(-20,0,0), 10, Vector(0.3, 0.5, 0.8), false, true);
    Object* S_transparent = new Sphere(Vector(0,0,0), 10, Vector(0.5, 0.5, 0.5), true, false);
    Object* S_transparent_hollow_outer = new Sphere(Vector(20,0,0), 10, Vector(0.7, 0.2, 0.5), true, false);
    // Object* S_transparent_hollow_inner = new Sphere(Vector(20,0,0), 9.8, Vector(0.7, 0.2, 0.5), true, false, true);
    Vector albedo(0.5, 0.5, 0.5); // grey sphere

    double I = 1E5; // intensity
    Vector light_position(-10, 20, 40);
    int nb_rays_per_pixel = 1;
    
    Object* S1 = new Sphere(Vector(0,0,1000), 940, Vector(0.9, 0.4, 0.3), false, false); // wall behind camera wall
    Object* S2 = new Sphere(Vector(0,-1000,0), 990, Vector(0.3, 0.4, 0.7), false, false); // floor
    Object* S3 = new Sphere(Vector(0,0,-1000), 940, Vector(0.4, 0.8, 0.7), false, false); // back wall
    Object* S4 = new Sphere(Vector(0,1000,0), 940, Vector(0.2, 0.5, 0.9), false, false); // ceiling
    Object* S5 = new Sphere(Vector(-1000,0,0), 940, Vector(0.9, 0.2, 0.9), false, false); // left wall
    Object* S6 = new Sphere(Vector(1000,0,0), 940, Vector(0.6, 0.5, 0.1), false, false); // right wall
    std::vector<Object*> objects; // = scene.objects;

    // objects.push_back(S_transparent);
    // objects.push_back(S_mirror);
    // objects.push_back(S_transparent_hollow_outer);
    // objects.push_back(S_transparent_hollow_inner);

    objects.push_back(S1);
    objects.push_back(S2);
    objects.push_back(S3);
    objects.push_back(S4);
    objects.push_back(S5);
    objects.push_back(S6);
    
    // Cat
    TriangleMesh cat_mesh(0.6, Vector(0, -10, 0), Vector(1., 1., 1.));
    cat_mesh.readOBJ("cat.obj");

    objects.push_back(&cat_mesh);

    Scene scene(objects);
 
    std::vector<unsigned char> image(W * H * 3, 0);
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Vector pixelColor = Vector(0., 0., 0.);
            double x, y;
            for (int k=0; k<nb_rays_per_pixel; k++) {
                boxMuller(0.5, x, y);
                Vector rand_dir = Vector(camera_origin[0] + (j+x) + 0.5 - W/2, camera_origin[1] - (i+y) - 0.5 + H/2, camera_origin[2] - W/(2*tan(fov/2)));  // as before but targetting pixel (i, j) + boxMuller() * spread
                rand_dir = rand_dir - camera_origin;
                rand_dir.normalize();
                Ray ray(camera_origin, rand_dir);
                pixelColor = pixelColor + scene.getColor(ray, nb_rays_per_pixel, scene, light_position, I);
            }
            Vector color = pixelColor;
            image[(i * W + j) * 3 + 0] = std::min(255., std::pow(color[0] /4., 1.0/2.2) * 255.);
            image[(i * W + j) * 3 + 1] = std::min(255., std::pow(color[1] /4., 1.0/2.2) * 255.);
            image[(i * W + j) * 3 + 2] = std::min(255., std::pow(color[2] /4., 1.0/2.2) * 255.);
        }
    }
    stbi_write_png("image.png", W, H, 3, &image[0], 0);

    // End the timer
    auto b = high_resolution_clock::now();
    std::cout << "Took " << duration_cast<seconds>(b - a).count() << " seconds" <<  std::endl;
    return 0;
}