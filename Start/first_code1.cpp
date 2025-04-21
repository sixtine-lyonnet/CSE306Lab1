#include <algorithm>
#include <chrono>
using namespace std::chrono;
#include <limits>
#include <iostream>
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

    Vector T1 = N;
    double min_value_index_N = 0; // std::min(abs(N[0]), abs(N[1]), abs(N[2]));
    for (int i = 1; i < 3; i++) {
        if (abs(N[i]) < abs(N[min_value_index_N])) {
            min_value_index_N = i;
        }
    }
    T1[min_value_index_N] = 0.;
    int index_swap_1 = abs(min_value_index_N - 1);
    int index_swap_2 = abs(min_value_index_N - 2);
    T1[index_swap_1] = N[index_swap_2];
    T1[index_swap_2] = -1 * N[index_swap_1];
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

class Intersection {
    public:
        Intersection(bool intersected, Vector P, Vector N, double t, int index) : intersected(intersected), P(P), N(N), t(t), index(index) {};

        bool intersected;
        Vector P;
        Vector N;
        double t;
        int index;
};

class Object {
    public:
        Object(const Vector& albedo = Vector(1., 1., 1.), bool isMirror = false, bool isTransparent = false): albedo(albedo), isMirror(isMirror), isTransparent(isTransparent) {}
        virtual Intersection intersect(const Ray& r) const = 0;

        Vector albedo;
        bool isMirror;
        bool isTransparent;
};

// _________________________________________________________________________________________________________________________________________________
// This code was taken from: https://pastebin.com/CAgp9r15

#include <string>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <vector>
 
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
    double scaling_factor;
    Vector translation;
    // Vector color;
    // double refractive_index;
    // bool reflects;
public:
    ~TriangleMesh() {}
    // TriangleMesh(double scaling_factor, Vector translation, Vector color = Vector(1., 1., 1.), double refractive_index = 1., bool reflects = false) : scaling_factor(scaling_factor), translation(translation), color(color), refractive_index(refractive_index), reflects(reflects) {};
    /*
    TriangleMesh(double scaling_factor, Vector translation, Vector albedo, double isMirror = false, bool isTransparent = false) : scaling_factor(scaling_factor), translation(translation) {
        this->albedo = albedo;
        this->isMirror = isMirror;
        this->isTransparent = isTransparent;
    };
    */
    TriangleMesh(double scaling_factor, Vector translation, Vector albedo, double isMirror = false, bool isTransparent = false) : Object(albedo, isMirror, isTransparent), scaling_factor(scaling_factor), translation(translation) {};


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
 
    }

    Intersection intersect(const Ray& r) const override {
    // bool intersect(const Ray& r, Vector &P, Vector &N, double &t) {
        Vector P;
        Vector N;
        double t = std::numeric_limits<double>::max();
        bool has_intersection = false;

        for (int i = 0; i < indices.size(); i++) {
            const Vector& A = scaling_factor * vertices[indices[i].vtxi] + translation;
            const Vector& B = scaling_factor * vertices[indices[i].vtxj] + translation;
            const Vector& C = scaling_factor * vertices[indices[i].vtxk] + translation;

            const Vector e1 = B - A;
            const Vector e2 = C - A;
            const Vector N = cross (e1, e2);
            // double invUN = 1. / dot(r.u, N);
            Vector AOcrossU = cross(A-r.origin, r.u);

            double localt = dot(A-r.origin, N) / dot(r.u, N); // implement the formula from the slides
            if (localt < 0) continue;
            if (localt > t) continue;

            double beta = dot(e2, AOcrossU) / dot(r.u, N); // implement formula from slides
            // check condition < 0 || > 1

            double gamma = dot(e1, AOcrossU) / dot(r.u, N); // implement the formula from the slides
            //check condition <0 || > 1

            double alpha = 1 - beta - gamma; // implement the formula from the slides
            // check condition <0 

            if (alpha >= 0.0f && alpha <= 1.0f && beta >= 0.0f && beta <= 1.0f && gamma >= 0.0f && gamma <= 1.0f && localt > 0.0f && localt < t) {
                t = localt;
                P = r.origin + t * r.u;
                Intersection intersection(true, P, N, t, -1);
                return intersection;
            }    
        }
        Intersection intersection(false, P, N, t, -1);
        return intersection;
    }
 
    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;
    
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

    Intersection intersect(const Ray& r) const override { // here, returned by reference but can do something else
        Vector P;
        Vector N;
        double t = std::numeric_limits<double>::max();
        double delta = sqr(dot(r.u, r.origin-C)) - ((r.origin-C).norm2() - sqr(R)); // discriminant formula from slides
        if (delta < 0) {
            Intersection intersection(false, P, N, t, -1);
            return intersection; // we only wants points in front of us?
        }
        double x = dot(r.u, C-r.origin);
        double t1 = x - sqrt(delta);
        double t2 = x + sqrt(delta);
        if (t2 < 0) {
            Intersection intersection(false, P, N, t, -1);
            return intersection; // there is not intersection
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
        
        Intersection intersection(true, P, N, t, -1);
        return intersection;
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

        Intersection intersect(const Ray& r) {
            Vector P;
            Vector N;
            double t;
            int id;
            bool result = false;
            double closest_t = std::numeric_limits<double>::max(); // maximum double so any sphere is closer initially
            
            for (int i = 0; i < objects.size(); i++) {

                // Vector localN, localP; // values for current local object
                // double localt;

                Intersection local_intersection = objects[i]->intersect(r);
                if (local_intersection.intersected) {
                    result = true;
                    if (local_intersection.t < closest_t) {
                        N = local_intersection.N;
                        P = local_intersection.P;
                        closest_t = local_intersection.t;
                        id = local_intersection.index;
                    }
                }
            }
            t = closest_t;
            Intersection intersection(result, P, N, t, id);
            return intersection;
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
        
            Intersection scene_intersection = scene.intersect(shadowRay);
            if (scene_intersection.intersected) {
                double intersectDist = (shadowP - scene_intersection.P).norm();
                if (intersectDist < lightDist) {
                    return true;
                }
            }
            return false;
        }


        Vector getColor(const Ray &ray, int bounce_number, Scene& scene, const Vector& light_position, const double& I) {
            if (bounce_number <= 0) {return Vector(0, 0, 0);} // ? // stop the iteration //else, decrease bounce_number
            // put nearly everything in it

            //Vector P, N;
            //double t; // maximum double so any sphere is closer initially;
            //int id;
            Vector color;

            Intersection scene_intersection = scene.intersect(ray);
            if (scene_intersection.intersected) {
                if (objects[scene_intersection.index]->isMirror) {
                    // Ray reflected_ray = ...;
                    Vector reflectDir = ray.u - 2 * dot(ray.u, scene_intersection.N) * scene_intersection.N;
                    Ray reflectedRay(scene_intersection.P + 0.00001 * scene_intersection.N, reflectDir);
                    return getColor(reflectedRay, bounce_number-1, scene, light_position, I);
                }

                if (objects[scene_intersection.index]->isTransparent) {
                    // Transparence
                    double n1 = 1;
                    double n2 = 1.5;
                    if (dot(ray.u, scene_intersection.N) > 0) {
                        std::swap(n1, n2);
                        scene_intersection.N = - scene_intersection.N;
                    }
                    /*
                    if (objects[id]->is_hollow) {
                        std::swap(n1, n2);
                        N = -N;
                    }
                    */
                    Vector lightDir = light_position - scene_intersection.P;
                    double lightDist = lightDir.norm();
                    Vector albedo = objects[scene_intersection.index]->albedo;
                    Vector omega_i = ray.u; // incoming ray;
                    Vector tangentialComponent = (n1 / n2) * (omega_i - dot(omega_i, scene_intersection.N) * scene_intersection.N); // formula from slides
                    double D = 1. - std::pow(n1/n2, 2)*(1. - std::pow(dot(omega_i, scene_intersection.N), 2)); // formula from slides // make sure stuff in 
                    Vector normalComponent;
                    if (D <= 0. ) {
                        // consider it to be a mirror
                        // normalComponent = -sqrt(-D) * N;
                        Vector reflectDir = ray.u - (2 * dot(ray.u, scene_intersection.N) * scene_intersection.N);
                        Ray reflectedRay(scene_intersection.P + 1e-4 * scene_intersection.N, reflectDir);
                        return getColor(reflectedRay, bounce_number-1, scene, light_position, I);
                    }
                    else {
                        normalComponent = -sqrt(D)*scene_intersection.N;
                        Vector transmittedDirection = normalComponent + tangentialComponent; //
                        transmittedDirection.normalize();

                        // Fresnel's law
                        double k0 = pow(n1 - n2, 2) / pow(n1 + n2, 2);
                        double R = k0 + (1 - k0) * pow(1 - abs(dot(scene_intersection.N, omega_i)), 5);
                        double T = 1- R;

                        // Solution (Leal gave me the hint)
                        Vector reflectDir = ray.u - (2 * dot(ray.u, scene_intersection.N) * scene_intersection.N);
                        Ray reflectedRay(scene_intersection.P + 1e-4 * scene_intersection.N, reflectDir);

                        Ray refracted_ray = Ray(scene_intersection.P - 1e-4 * scene_intersection.N, transmittedDirection);

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
                    
                    if (isShadow(scene, light_position, scene_intersection.P, scene_intersection.N, scene_intersection.t)) {
                        visibility = 0.;
                    }
                    
                    Vector lightDir = light_position - scene_intersection.P;
                    double lightDist = lightDir.norm();
                    Vector albedo = objects[scene_intersection.index]->albedo;
                    Vector omega_i = lightDir/lightDir.norm();
                    double dot_product = std::max(dot(scene_intersection.N, omega_i), 0.); // make sure not negative
                    Lo = (I/(4 * PI * lightDir.norm2())) * (albedo / PI) * visibility * dot_product;

                    // add indirect lighting
                    // Ray randomRay = Ray(P, random_cos(N));; // randomly sample ray using random_cos
                    // Lo = Lo + albedo * getColor(randomRay, bounce_number-1, scene, light_position, I);
                    
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
    int W = 64; //512;
    int H = 64; // 512;
    Vector camera_origin(0, 0, 55);
    double fov = 60 * PI / 180; // make sure radians not degrees
    /*
    Object* S_mirror = new Sphere(Vector(-20,0,0), 10, Vector(0.3, 0.5, 0.8), false, true);
    Object* S_transparent = new Sphere(Vector(0,0,0), 10, Vector(0.5, 0.5, 0.5), true, false);
    Object* S_transparent_hollow_outer = new Sphere(Vector(20,0,0), 10, Vector(0.7, 0.2, 0.5), true, false);
    Object* S_transparent_hollow_inner = new Sphere(Vector(20,0,0), 9.8, Vector(0.7, 0.2, 0.5), true, false, true);
    Vector albedo(0.5, 0.5, 0.5); // grey sphere
    */
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
    /*
    objects.push_back(S_transparent);
    objects.push_back(S_mirror);
    objects.push_back(S_transparent_hollow_outer);
    objects.push_back(S_transparent_hollow_inner);
    */
    objects.push_back(S1);
    objects.push_back(S2);
    objects.push_back(S3);
    objects.push_back(S4);
    objects.push_back(S5);
    objects.push_back(S6);
    
    // Cat
    // TriangleMesh cat_mesh(0.6, Vector(0, -10, 0), Vector(1., 1., 1.));
    // cat_mesh.readOBJ("cat.obj");

    //TriangleMesh* cat_mesh = new TriangleMesh(0.6, Vector(0, -10, 0), Vector(1., 1., 1.));
    TriangleMesh cat_mesh(0.6, Vector(0, -10, 0), Vector(1., 1., 1.));
    //cat_mesh->readOBJ("cat.obj");
    cat_mesh.readOBJ("cat.obj");

    //objects.push_back(cat_mesh);
    objects.push_back(&cat_mesh);

    Scene scene(objects);
 
    std::vector<unsigned char> image(W * H * 3, 0);
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            // double d = -W/(2*tan(fov/2));
            // Vector ray_direction(j-W/2+0.5, H/2-i+0.5, d);
            // ray_direction.normalize();
            // Ray r(camera_origin, ray_direction);
            Vector pixelColor = Vector(0., 0., 0.);
            double x, y;
            for (int k=0; k<nb_rays_per_pixel; k++) {
                boxMuller(0.5, x, y);
                // double d = -W/(2*tan(fov/2));
                // Vector ray_direction(j-W/2+0.5, H/2-i+0.5, d);
                Vector rand_dir = Vector(camera_origin[0] + (j+x) + 0.5 - W/2, camera_origin[1] - (i+y) - 0.5 + H/2, camera_origin[2] - W/(2*tan(fov/2)));  // as before but targetting pixel (i, j) + boxMuller() * spread
                rand_dir = rand_dir - camera_origin;
                rand_dir.normalize();
                Ray ray(camera_origin, rand_dir);
                pixelColor = pixelColor + scene.getColor(ray, nb_rays_per_pixel, scene, light_position, I);
            }
            Vector color = pixelColor;
            // Vector color = scene.getColor(r, 2, scene, light_position, I);
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