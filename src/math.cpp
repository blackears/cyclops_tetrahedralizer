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


#include "math.h"

using namespace CyclopsTetra3D;

const Vector3 Vector3::X_POS = Vector3(1, 0, 0);
const Vector3 Vector3::X_NEG = Vector3(-1, 0, 0);
const Vector3 Vector3::Y_POS = Vector3(0, 1, 0);
const Vector3 Vector3::Y_NEG = Vector3(0, -1, 0);
const Vector3 Vector3::Z_POS = Vector3(0, 0, 1);
const Vector3 Vector3::Z_NEG = Vector3(0, 0, -1);

//Vector3 tetrahedron_circumcenter(const Vector3& p0, const Vector3& p1, const Vector3& p2, const Vector3& p3) {
//    //https://rodolphe-vaillant.fr/entry/127/find-a-tetrahedron-circumcenter
//
//    //From Matthias Muller
//    //https://github.com/matthias-research/pages/blob/62fa5a972572338a9afb7f50bfd22aa8d7d90e19/tenMinutePhysics/BlenderTetPlugin.py#L68
//    Vector3 b = p1 - p0;
//    Vector3 c = p2 - p0;
//    Vector3 d = p3 - p0;
//
//    real det = 2.0 * (b.x * (c.y * d.z - c.z * d.y)
//        - b.y * (c.x * d.z - c.z * d.x)
//        + b.z * (c.x * d.y - c.y * d.x));
//
//    if (det == 0.0) {
//        return p0;
//    }
//    else {
//        Vector3 v = c.cross(d) * b.dot(b) + d.cross(b) * c.dot(c) + b.cross(c) * d.dot(d);
//        v /= det;
//        return p0 + v;
//    }
//}
