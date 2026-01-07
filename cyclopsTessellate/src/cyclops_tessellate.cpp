/*
 * Copyright (c) 2026 Mark McKay
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */


#include "cyclops_tessellate.h"

using namespace std;
using namespace CyclopsTessellate3D;

//https://github.com/matthias-research/pages/blob/62fa5a972572338a9afb7f50bfd22aa8d7d90e19/tenMinutePhysics/BlenderTetPlugin.py

template<typename T>
Vector3<T> Tetrahedron<T>::calc_circum_center(Vector3<T> p0, Vector3<T> p1, Vector3<T> p2, Vector3<T> p3) const {
    //https://rodolphe-vaillant.fr/entry/127/find-a-tetrahedron-circumcenter

    //From Matthias Muller
    //https://github.com/matthias-research/pages/blob/62fa5a972572338a9afb7f50bfd22aa8d7d90e19/tenMinutePhysics/BlenderTetPlugin.py#L68

    Vector3<T> b = p1 - p0;
    Vector3<T> c = p2 - p0;
    Vector3<T> d = p3 - p0;
 
    T det = 2.0 * (b.x * (c.y * d.z - c.z * d.y) 
        - b.y * (c.x * d.z - c.z * d.x) 
        + b.z * (c.x * d.y - c.y * d.x));

    if (det == 0.0) {
        return p0;
    }
    else {
        Vector3<T> v = c.cross(d) * b.dot(b) + d.cross(b) * c.dot(c) + b.cross(c) * d.dot(d);
        v /= det;
        return p0 + v;
    }
}

//return value on [0 - 1] where 1 is a perfect tetrahedron
template<typename T>
T Tetrahedron<T>::quality(Vector3<T> p0, Vector3<T> p1, Vector3<T> p2, Vector3<T> p3) const {
    Vector3<T> d0 = p1 - p0;
    Vector3<T> d1 = p2 - p0;
    Vector3<T> d2 = p3 - p0;
    Vector3<T> d3 = p2 - p1;
    Vector3<T> d4 = p3 - p2;
    Vector3<T> d5 = p1 - p3;

    T s0 = d0.magnitude();
    T s1 = d1.magnitude();
    T s2 = d2.magnitude();
    T s3 = d3.magnitude();
    T s4 = d4.magnitude();
    T s5 = d5.magnitude();

    T ms = (s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3 + s4 * s4 + s5 * s5) / 6.0;
    T rms = sqrt(ms);

    T s = 12.0 / sqrt(2.0);

    T vol = d0.dot(d1.cross(d2)) / 6.0;
    return s * vol / (rms * rms * rms);
}


// Vector3<T> Tetrahedron<T>::calc_circum_center(Vector3<T> v0, Vector3<T> v1, Vector3<T> v2, Vector3<T> v3) const {
//     //https://rodolphe-vaillant.fr/entry/127/find-a-tetrahedron-circumcenter

//     //Row vectors
//     Matrix3x3 A(v1 - v0,
//         v2 - v0,
//         v3 - v0);
    
//     Vector3 b(v1.dot(v1) - v0.dot(v0), 
//         v2.dot(v2) - v0.dot(v0),
//         v3.dot(v3) - v0.dot(v0)) / 2.0;

//     Vector c = A.inverse() * b;
    
//     return c;
// }

// template<typename T>
// bool BVHTree<T>::is_inside(Vector3<T> p, T dist_min) const {
//     return true;
// }

template<typename T>
void BVHTree<T>::build_from_triangles(Vector3<T>* tri_points) {

}

template<typename T>
bool BVHTree<T>::is_inside(Vector3<T> p, T dist_min) const {
    //Check multiple directions to minimize chance of hitting an edge
    const Vector3[T] dirs = {
        Vector3<T>(1,0,0),
        Vector3<T>(-1,0,0),
        Vector3<T>(0,1,0),
        Vector3<T>(0,-1,0),
        Vector3<T>(0,0,1),
        Vector3<T>(0,0,-1)
    };

    int num_inside = 0;

    for (int i = 0; i < 6; i++) {
        Vector3<T> hit_pos;
        Vector3<T> hit_normal;
        int hit_index;
        T hit_distance;

        this->ray_cast(p, dirs[i], 1e10, hit_pos, hit_normal, hit_index, hit_distance);
        
        //Check hit is valid
        if (hit_index >= 0) {
            if (hit_normal.dot(dirs[i]) > 0) {
                num_inside++;
            }

            if (dist_min > 0.0 && hit_distance < dist_min) {
                //If hit outside face, we are outside
                return false;
            }
        }
    }

    return num_inside >= 3;
}

template<typename T>
void BVHTree<T>::ray_cast(Vector3<T> origin, Vector3<T> direction, T distance, Vector3<T> &hit_pos, Vector3<T> &hit_normal, int &index, T &hit_distance) const {

}


template<typename T>
void CyclopsTess3D<T>::tessellate_tetrahedra(float* points) {
}