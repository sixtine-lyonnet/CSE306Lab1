#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
 
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
 
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define PI 3.14159

double sqr(double x) {return x * x;};
 
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
Vector operator*(const double a, const Vector& b) {
    return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const double b) {
    return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator/(const Vector& a, const double b) {
    return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

class Ray {
public:
    Ray(const Vector& O, const Vector &u): origin(O), u(u) {};
    Vector origin;
    Vector u;
    //void normalize() {
    //    u.normalize();
    //}

};
 
class Sphere {
public:
    Sphere(const Vector& C, double R) : C(C), R(R) {};
    Vector C;
    double R;

    bool intersect(const Ray& r) {
        double delta = sqr(dot(r.u, r.origin-C)) - (r.origin-C).norm2() + sqr(R); // discriminant formula from slides
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
        return true;
    }
};
 
int main() {
    int W = 512;
    int H = 512;
    Vector camera_origin(0, 0, 55);
    double fov = 60 * PI / 180; // make sure radians not degrees
    Sphere S(Vector(0,0,0), 10);
 
    std::vector<unsigned char> image(W * H * 3, 0);
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {

            double d = -W/(2*tan(fov/2));
            Vector ray_direction(j-W/2+0.5, H/2-i+0.5, d);
            ray_direction.normalize();
            Ray r(camera_origin, ray_direction);
            if (S.intersect(r)) {
                image[(i * W + j) * 3 + 0] = 255;
                image[(i * W + j) * 3 + 1] = 255;
                image[(i * W + j) * 3 + 2] = 255;
            }
 
            
        }
    }
    stbi_write_png("image.png", W, H, 3, &image[0], 0);
 
    return 0;
}