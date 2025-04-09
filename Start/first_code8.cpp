#include <algorithm>
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
    T1[index_swap_2] = - 1 * N[index_swap_1]; // minus
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
 
class Sphere {
public:
    Sphere(const Vector& C, double R, Vector albedo, bool is_T, bool is_M, bool is_H) : C(C), R(R), albedo(albedo), transparent(is_T), mirror(is_M), is_hollow(is_H) {};
    Sphere(const Vector& C, double R, Vector albedo, bool is_T, bool is_M) : C(C), R(R), albedo(albedo), transparent(is_T), mirror(is_M), is_hollow(false) {};
    Vector C;
    double R;
    // double t;
    Vector albedo;
    bool transparent;
    bool mirror;
    bool is_hollow;

    bool intersect(const Ray& r, Vector &P, Vector &N, double &t) { // here, returned by reference but can do something else
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

class Intersection {
    public:
        Intersection(bool intersected, Vector P, Vector N, double t, int index) : intersected(intersected), P(P), N(N), t(t), index(index) {};

        bool intersected;
        Vector P;
        Vector N;
        double t;
        int index;
};

class Scene {
    public:
        Scene(std::vector<Sphere> obj) : objects(obj) {};
        
        std::vector<Sphere> objects;
        // add albedo
        // add light color?

        void add(const Sphere &s) {
            objects.push_back(s);
        }

        bool intersect(const Ray& r, Vector& P, Vector& N, double& t, int &id) {
            bool result = false;
            double closest_t = std::numeric_limits<double>::max(); // maximum double so any sphere is closer initially
            
            for (int i = 0; i < objects.size(); i++) {

                Vector localN, localP; // values for current local object
                double localt;

                if (objects[i].intersect(r, localP, localN, localt)) {
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
                if (objects[id].mirror) {
                    // Ray reflected_ray = ...;
                    Vector reflectDir = ray.u - 2 * dot(ray.u, N) * N;
                    Ray reflectedRay(P + 0.00001 * N, reflectDir);
                    return getColor(reflectedRay, bounce_number-1, scene, light_position, I);
                }

                if (objects[id].transparent) {
                    // Transparence
                    double n1 = 1;
                    double n2 = 1.5;
                    if (dot(ray.u, N) > 0) {
                        std::swap(n1, n2);
                        N = -N;
                    }
                    if (objects[id].is_hollow) {
                        std::swap(n1, n2);
                        N = -N;
                    }
                    Vector lightDir = light_position - P;
                    double lightDist = lightDir.norm();
                    Vector albedo = objects[id].albedo;
                    Vector omega_i = ray.u; // incoming ray;
                    Vector tangentialComponent = (n1 / n2) * (omega_i - dot(omega_i, N) * N); // formula from slides
                    double D = 1 - std::pow(n1/n2, 2)*(1-std::pow(dot(omega_i, N), 2)); // formula from slides // make sure stuff in 
                    Vector normalComponent;
                    if (D < 0) {
                        // consider it to be a mirror (copy 3 lines)
                        // normalComponent = -sqrt(-D) * N;
                        Vector reflectDir = ray.u - 2 * dot(ray.u, N) * N;
                        Ray reflectedRay(P + 0.00001 * N, reflectDir);
                        return getColor(reflectedRay, bounce_number-1, scene, light_position, I);
                    }
                    else {
                        normalComponent = -sqrt(D)*N;
                        Vector transmittedDirection = normalComponent + tangentialComponent; //
                        transmittedDirection.normalize();
                        Ray transmittedRay(P - 1e-4 * N, transmittedDirection); // go inside object (make sure you don't interest with same border)
                        Vector transmittedColor = getColor(transmittedRay, bounce_number-1, scene, light_position, I); //
                        return transmittedColor;
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
                    Vector albedo = objects[id].albedo;
                    Vector omega_i = lightDir/lightDir.norm();
                    double dot_product = std::max(dot(N, omega_i), 0.); // make sure not negative
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
    /*
        Vector del1 = Vector(2, -2, 1);
        std::cout << "del1 = " << del1[0] << ", " << del1[1] << ", " << del1[2] << std::endl;
        Vector del2 = random_cos(del1);
        std::cout << "del2 = " << del2[0] << ", " << del2[1] << ", " << del2[2] << std::endl;
    */
    int W = 512;
    int H = 512;
    Vector camera_origin(0, 0, 55);
    double fov = 60 * PI / 180; // make sure radians not degrees
    Sphere S_transparent(Vector(0,0,0), 10, Vector(0.5, 0.5, 0.5), true, false);
    Sphere S_mirror(Vector(-20,0,0), 10, Vector(0.3, 0.5, 0.8), false, true);
    Sphere S_transparent_hollow_outer(Vector(20,0,0), 10, Vector(0.7, 0.2, 0.5), true, false);
    Sphere S_transparent_hollow_inner(Vector(20,0,0), 9.8, Vector(0.7, 0.2, 0.5), true, false, true);
    Vector albedo(0.5, 0.5, 0.5); // grey sphere
    double I = 1E5; // intensity
    Vector light_position(-10, 20, 40);
    
    Sphere S1(Vector(0,0,1000), 940, Vector(0.9, 0.4, 0.3), false, false); // wall behind camera wall
    Sphere S2(Vector(0,-1000,0), 990, Vector(0.3, 0.4, 0.7), false, false); // floor
    Sphere S3(Vector(0,0,-1000), 940, Vector(0.4, 0.8, 0.7), false, false); // back wall
    Sphere S4(Vector(0,1000,0), 940, Vector(0.2, 0.5, 0.9), false, false); // ceiling
    Sphere S5(Vector(-1000,0,0), 940, Vector(0.9, 0.2, 0.9), false, false); // left wall
    Sphere S6(Vector(1000,0,0), 940, Vector(0.6, 0.5, 0.1), false, false); // right wall
    std::vector<Sphere> objects; // = scene.objects;
    objects.push_back(S_transparent);
    objects.push_back(S_mirror);
    objects.push_back(S_transparent_hollow_outer);
    objects.push_back(S_transparent_hollow_inner);
    objects.push_back(S1);
    objects.push_back(S2);
    objects.push_back(S3);
    objects.push_back(S4);
    objects.push_back(S5);
    objects.push_back(S6);
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
            for (int k=0; k<4; k++) {
                boxMuller(0.5, x, y);
                // double d = -W/(2*tan(fov/2));
                // Vector ray_direction(j-W/2+0.5, H/2-i+0.5, d);
                Vector rand_dir = Vector(camera_origin[0] + (j+x) + 0.5 - W/2, camera_origin[1] - (i+y) - 0.5 + H/2, camera_origin[2] - W/(2*tan(fov/2)));  // as before but targetting pixel (i, j) + boxMuller() * spread
                rand_dir = rand_dir - camera_origin;
                rand_dir.normalize();
                Ray ray(camera_origin, rand_dir);
                pixelColor = pixelColor + scene.getColor(ray, 4, scene, light_position, I);
            }
            Vector color = pixelColor;
            // Vector color = scene.getColor(r, 4, scene, light_position, I);
            image[(i * W + j) * 3 + 0] = std::min(255., std::pow(color[0] /4., 1.0/2.2) * 255.);
            image[(i * W + j) * 3 + 1] = std::min(255., std::pow(color[1] /4., 1.0/2.2) * 255.);
            image[(i * W + j) * 3 + 2] = std::min(255., std::pow(color[2] /4., 1.0/2.2) * 255.);
        }
    }
    stbi_write_png("image8.png", W, H, 3, &image[0], 0);
 
    return 0;
}