
#ifndef SYSTEMSTRUCTURES_H_
#define SYSTEMSTRUCTURES_H_

#include <memory>
#include <cmath>
#include <random>
#include <algorithm>
#include <iterator>
#include <thrust/transform_reduce.h>
#include <thrust/functional.h>

#include <thrust/random.h>
#include <thrust/extrema.h>
#include <thrust/sort.h>
#include <thrust/for_each.h>
#include <thrust/transform.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/tuple.h>
#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
#include <thrust/remove.h>
#include <thrust/copy.h>
#include <thrust/execution_policy.h>
#include <thrust/pair.h>
#include <thrust/unique.h>
#include <thrust/remove.h>
#include <thrust/binary_search.h>
#include <thrust/reduce.h>
#include <thrust/replace.h>
#include <thrust/gather.h>
#include <thrust/random/uniform_real_distribution.h>
#include <stdint.h>
#include <thrust/sequence.h>
#include <thrust/host_vector.h>

// Type Aliases
typedef thrust::tuple<int, bool, double> Tubd;
typedef thrust::tuple<int, bool> Tub;
typedef thrust::tuple<int, double> Tud;
typedef thrust::tuple<bool, double> Tbd;

typedef thrust::tuple<int, int, double> Tuud;

typedef thrust::tuple<int, int, int, int, int, int, int> Tuuuuuuu;
typedef thrust::tuple<int, int, int, int, double> Tuuuud;
typedef thrust::tuple<int, int, int, int, int> Tuuuuu;
typedef thrust::tuple<int, int, int, int> Tuuuu;
typedef thrust::tuple<int, int, int, double> Tuuud;
typedef thrust::tuple<int, int, int> Tuuu;
typedef thrust::tuple<int, int> Tuu;

typedef thrust::tuple<int, double, double, double> Tuddd;
typedef thrust::tuple<double, double, double, int> Tdddu;
typedef thrust::tuple<double, double> Tdd;

typedef thrust::tuple<bool, double, double, double, double, double, double> BoolCVec6;
typedef thrust::tuple<int, double, double, double, double, double, double,double, double, double> UCVec9;

typedef thrust::tuple<int, double, double, double, double, double, double> UCVec6;
typedef thrust::tuple<int, int, double, double, double, double> U2CVec4;
typedef thrust::tuple<int, int, double, double, double> U2CVec3;
typedef thrust::tuple<int, double, double, double> UCVec3;
typedef thrust::tuple<bool, double, double, double> BoolCVec3;

typedef thrust::tuple<double, double, double, double, double, double, double, double, double> CVec9;
typedef thrust::tuple<double, double, double, double, double, double, double, double> CVec8;
typedef thrust::tuple<double, double, double, double, double, double, double> CVec7;
typedef thrust::tuple<double, double, double, double, double, double> CVec6;
typedef thrust::tuple<double, double, double, double, double> CVec5;
typedef thrust::tuple<double, double, double, double> CVec4;
typedef thrust::tuple<double, double, double> CVec3;
typedef thrust::tuple<double, double> CVec2;

typedef thrust::tuple<CVec3, CVec3, CVec3> Mat_3x3;

// Forward Declarations of Structs
struct GeneralParams;
struct DomainParams;
struct AuxVecs;
struct CoordInfoVecs;
struct CapsidInfoVecs;
struct BendingTriangleInfoVecs;
struct AreaTriangleInfoVecs;
struct LinearSpringInfoVecs;
struct LJInfoVecs;


struct HostSetInfoVecs {

    // Vectors to store information about nodes, triangles, edges in the upper hemisphere.
    thrust::host_vector<double> scaling_per_edge;
    thrust::host_vector<int> nodes_in_upperhem;
    thrust::host_vector<int> triangles_in_upperhem;
    thrust::host_vector<int> edges_in_upperhem;
    thrust::host_vector<int> edges_in_upperhem_list;
    thrust::host_vector<int> boundaries_in_upperhem;
    thrust::host_vector<int> boundaries_in_lowerhem;
    thrust::host_vector<int> nodes_in_tip;
    thrust::host_vector<int> edges_in_tip;
    thrust::host_vector<int> triangles_in_tip; // redundant
    thrust::host_vector<int> nodes_in_DV; // 0 for outDV regions, 1 for inDV regions
    thrust::host_vector<int> edges_in_DV; // 0 for outDV regions, 1 for inDV regions   

    // Vectors to store the relationship between nodes and triangles (local coordinates).
    thrust::host_vector<int> nodes2Triangles_1;
    thrust::host_vector<int> nodes2Triangles_2;
    thrust::host_vector<int> nodes2Triangles_3;
    thrust::host_vector<int> nodes2Triangles_4;
    thrust::host_vector<int> nodes2Triangles_5;
    thrust::host_vector<int> nodes2Triangles_6;
    thrust::host_vector<int> nodes2Triangles_7;
    thrust::host_vector<int> nodes2Triangles_8;
    thrust::host_vector<int> nodes2Triangles_9;

    // Vectors to store neighborhood data for nodes.
    thrust::host_vector<int> nndata1;
    thrust::host_vector<int> nndata2;
    thrust::host_vector<int> nndata3;
    thrust::host_vector<int> nndata4;
    thrust::host_vector<int> nndata5;
    thrust::host_vector<int> nndata6;
    thrust::host_vector<int> nndata7;
    thrust::host_vector<int> nndata8;
    thrust::host_vector<int> nndata9;

    // Vectors to store coordinates and forces for nodes.
    thrust::host_vector<double> capsidNodeLocX;
    thrust::host_vector<double> capsidNodeLocY;
    thrust::host_vector<double> capsidNodeLocZ;
    thrust::host_vector<bool> isNodeFixed;
    thrust::host_vector<double> nodeLocX;
    thrust::host_vector<double> nodeLocY;
    thrust::host_vector<double> nodeLocZ;
    thrust::host_vector<double> nodeForceX;
    thrust::host_vector<double> nodeForceY;
    thrust::host_vector<double> nodeForceZ;
    thrust::host_vector<double> nodeVelX;
    thrust::host_vector<double> nodeVelY;
    thrust::host_vector<double> nodeVelZ;

    // Vectors to store local coordinates of triangles.
    thrust::host_vector<int> triangles2Nodes_1;
    thrust::host_vector<int> triangles2Nodes_2;
    thrust::host_vector<int> triangles2Nodes_3;
    thrust::host_vector<int> triangles2Triangles_1;
    thrust::host_vector<int> triangles2Triangles_2;
    thrust::host_vector<int> triangles2Triangles_3;

    // Vectors to store local coordinates of edges.
    thrust::host_vector<int> edges2Nodes_1;
    thrust::host_vector<int> edges2Nodes_2;
    thrust::host_vector<int> edges2Triangles_1;
    thrust::host_vector<int> edges2Triangles_2;

    // Vectors to store local coordinates of edges on triangles.
    thrust::host_vector<int> triangles2Edges_1;
    thrust::host_vector<int> triangles2Edges_2;
    thrust::host_vector<int> triangles2Edges_3;

    // Vector to store the initial length of edges.
    thrust::host_vector<double> edge_initial_length;
    thrust::host_vector<double> edge_rest_length;
    
//    thrust::host_vector<double> pathlength_scaled;
//    
//    // Radial basis vector (in-surface direction pointing away from center)
//    thrust::host_vector<double> e_R_x;
//    thrust::host_vector<double> e_R_y;
//    thrust::host_vector<double> e_R_z;
//    
//    // Circumferential basis vector (in-surface, perpendicular to e_R)
//    thrust::host_vector<double> e_phi_x;
//    thrust::host_vector<double> e_phi_y;
//    thrust::host_vector<double> e_phi_z;
//    
//    // Height/normal basis vector (surface normal, points radially outward)
//    thrust::host_vector<double> e_h_x;
//    thrust::host_vector<double> e_h_y;
//    thrust::host_vector<double> e_h_z;
};

// Functor to add forces to nodes in the system.
struct AddForceFunctor {
    double* forceXAddr;
    double* forceYAddr;
    double* forceZAddr;

    // Constructor to set the addresses of force vectors.
    __host__ __device__
    AddForceFunctor(double* _forceXAddr, double* _forceYAddr, double* _forceZAddr) :
        forceXAddr(_forceXAddr), forceYAddr(_forceYAddr), forceZAddr(_forceZAddr) {}

    // Overloaded operator to add forces to nodes.
    __device__
    void operator() (const Tuddd& u1d3) {
        int idToAssign = thrust::get<0>(u1d3);
        if (!isnan(thrust::get<1>(u1d3)) && !isnan(thrust::get<2>(u1d3)) && !isnan(thrust::get<3>(u1d3))) {
            forceXAddr[idToAssign] += thrust::get<1>(u1d3);
            forceYAddr[idToAssign] += thrust::get<2>(u1d3);
            forceZAddr[idToAssign] += thrust::get<3>(u1d3);
        }
    }
};

// Functor to add forces to nodes in the system (with an additional check for the maximum node count).
struct AddForceFunctorAlt {
    int maxNodeCount;
    double* forceXAddr;
    double* forceYAddr;
    double* forceZAddr;

    // Constructor to set the maximum node count and addresses of force vectors.
    __host__ __device__
    AddForceFunctorAlt(int& _maxNodeCount, double* _forceXAddr, double* _forceYAddr, double* _forceZAddr) :
        maxNodeCount(_maxNodeCount), forceXAddr(_forceXAddr), forceYAddr(_forceYAddr), forceZAddr(_forceZAddr) {}

    // Overloaded operator to add forces to nodes.
    __device__
    void operator() (const Tuddd& u1d3) {
        int idToAssign = thrust::get<0>(u1d3);
        if (idToAssign < maxNodeCount) {
            forceXAddr[idToAssign] += thrust::get<1>(u1d3);
            forceYAddr[idToAssign] += thrust::get<2>(u1d3);
            forceZAddr[idToAssign] += thrust::get<3>(u1d3);
        }
    }
};

//#endif // SYSTEMSTRUCTURES_H_


// Functor for element-wise addition of two UCVec3 tuples. This version is deprecated. New function to be defined. 
//struct UCVec3Add : public thrust::binary_function<UCVec3, UCVec3, UCVec3> {
//    __host__ __device__
//    UCVec3 operator()(const UCVec3& vec1, const UCVec3& vec2) {
//        return thrust::make_tuple(
//            thrust::get<0>(vec1) + thrust::get<0>(vec2),
//            thrust::get<1>(vec1) + thrust::get<1>(vec2),
//            thrust::get<2>(vec1) + thrust::get<2>(vec2),
//            thrust::get<3>(vec1) + thrust::get<3>(vec2)
//        );
//    }
//};

struct UCVec3Add {
    using first_argument_type = UCVec3;
    using second_argument_type = UCVec3;
    using result_type = UCVec3;
    
    __host__ __device__ UCVec3 operator() (const UCVec3& vec1, const UCVec3& vec2){
        return thrust::make_tuple(
            thrust::get<0>(vec1) + thrust::get<0>(vec2),
            thrust::get<1>(vec1) + thrust::get<1>(vec2),
            thrust::get<2>(vec1) + thrust::get<2>(vec2),
            thrust::get<3>(vec1) + thrust::get<3>(vec2)
        );
    }
};

//// Functor for element-wise addition of two CVec3 tuples.
//struct CVec3Add : public thrust::binary_function<CVec3, CVec3, CVec3> {
//    __host__ __device__
//    CVec3 operator()(const CVec3& vec1, const CVec3& vec2) {
//        return thrust::make_tuple(
//            thrust::get<0>(vec1) + thrust::get<0>(vec2),
//            thrust::get<1>(vec1) + thrust::get<1>(vec2),
//            thrust::get<2>(vec1) + thrust::get<2>(vec2)
//        );
//    }
//};

struct CVec3Add {
    using first_argument_type = CVec3;
    using second_argument_type = CVec3;
    using result_type = CVec3;
    
    __host__ __device__ CVec3 operator()(const CVec3& vec1, const CVec3& vec2){
        return thrust::make_tuple(
            thrust::get<0>(vec1) + thrust::get<0>(vec2),
            thrust::get<1>(vec1) + thrust::get<1>(vec2),
            thrust::get<2>(vec1) + thrust::get<2>(vec2)            
        );
    }
};

//// Functor for element-wise addition of two CVec4 tuples.
//struct CVec4Add : public thrust::binary_function<CVec4, CVec4, CVec4> {
//    __host__ __device__
//    CVec4 operator()(const CVec4& vec1, const CVec4& vec2) {
//        return thrust::make_tuple(
//            thrust::get<0>(vec1) + thrust::get<0>(vec2),
//            thrust::get<1>(vec1) + thrust::get<1>(vec2),
//            thrust::get<2>(vec1) + thrust::get<2>(vec2),
//            thrust::get<3>(vec1) + thrust::get<3>(vec2)
//        );
//    }
//};

struct CVec4Add {
    using first_argument_type = CVec4;
    using second_argument_type = CVec4;
    using result_type = CVec4;
    
    __host__ __device__ CVec4 operator()(const CVec4& vec1, const CVec4& vec2){
         return thrust::make_tuple(
            thrust::get<0>(vec1) + thrust::get<0>(vec2),
            thrust::get<1>(vec1) + thrust::get<1>(vec2),
            thrust::get<2>(vec1) + thrust::get<2>(vec2),
            thrust::get<3>(vec1) + thrust::get<3>(vec2)
         );
    }
};

// Functor for calculating the norm (magnitude) of a CVec3 tuple as a binary function.
struct CVec3NormBinary {
    __host__ __device__
    double operator()(const CVec3& vec1, const CVec3& vec2) {
        // Calculate the Euclidean distance between vec1 and vec2.
        // The magnitude of the resulting CVec3 is returned as the norm.
        return fabs(
            ((thrust::get<0>(vec1) - thrust::get<0>(vec2))) +
            ((thrust::get<1>(vec1) - thrust::get<1>(vec2))) +
            ((thrust::get<2>(vec1) - thrust::get<2>(vec2)))
        );
    }
};

// Functor for calculating the norm (magnitude) of a CVec3 tuple as a unary function.
struct CVec3NormUnary {
    __host__ __device__
    double operator()(const CVec3& vec) {
        // Calculate the Euclidean distance of the CVec3 tuple from the origin (0,0,0).
        // The magnitude of the CVec3 is returned as the norm.
        return (
            sqrt(thrust::get<0>(vec) * thrust::get<0>(vec) +
            thrust::get<1>(vec) * thrust::get<1>(vec) +
            thrust::get<2>(vec) * thrust::get<2>(vec))
        );
    }
};

///////////////////////////////////////////////
// Random number generators

// Functor for generating normally distributed random numbers.
struct psrnormgen {
    double a, b;

    __host__ __device__
    psrnormgen(double _a, double _b) : a(_a), b(_b) {}

    __device__ double operator()(const int n) const {
        thrust::default_random_engine rng(n);
        thrust::normal_distribution<double> dist(a, b);
        rng.discard(n);
        return dist(rng);
    }
};

// Functor for generating uniformly distributed random numbers.
struct psrunifgen {
    double a, b;

    __host__ __device__
    psrunifgen(double _a, double _b) : a(_a), b(_b) {}

    __device__ double operator()(const int n) const {
        thrust::default_random_engine rng(n);
        thrust::uniform_real_distribution<double> dist(a, b);
        rng.discard(n);
        return dist(rng);
    }
};

// Functor for calculating the torsion angle of a triplet of nodes.
struct TorsionAngleFunctor {
    double* locXAddr;
    double* locYAddr;
    double* locZAddr;

    __host__ __device__
    TorsionAngleFunctor(double* _locXAddr, double* _locYAddr, double* _locZAddr) :
        locXAddr(_locXAddr), locYAddr(_locYAddr), locZAddr(_locZAddr) {}

    __device__
    double operator()(const Tuuu& u3) {
        int indexLeft = thrust::get<0>(u3);
        int indexCenter = thrust::get<1>(u3);
        int indexRight = thrust::get<2>(u3);

        // Calculate the distance between left and center nodes.
        double distLCX = locXAddr[indexLeft] - locXAddr[indexCenter];
        double distLCY = locYAddr[indexLeft] - locYAddr[indexCenter];
        double distLCZ = locZAddr[indexLeft] - locZAddr[indexCenter];

        // Calculate the distance between right and center nodes.
        double distRCX = locXAddr[indexRight] - locXAddr[indexCenter];
        double distRCY = locYAddr[indexRight] - locYAddr[indexCenter];
        double distRCZ = locZAddr[indexRight] - locZAddr[indexCenter];

        // Calculate the lengths between left & center and right & center.
        double lenLC = sqrt(distLCX * distLCX + distLCY * distLCY + distLCZ * distLCZ);
        double lenRC = sqrt(distRCX * distRCX + distRCY * distRCY + distRCZ * distRCZ);

        // Calculate the normalized dot product to obtain the cosine of the angle.
        double cosTheta = (distLCX * distRCX + distLCY * distRCY + distLCZ * distRCZ) / (lenLC * lenRC);

        // Handling rounding errors for acos.
        if (cosTheta < -1.0)
            cosTheta = -1.0;
        else if (cosTheta > 1.0)
            cosTheta = 1.0;

        // Calculate the torsion angle (radians) using the inverse cosine (acos).
        return acos(cosTheta);
    }
};

// Functor for averaging strain values.
struct AveStrainFunctor {
    __host__ __device__
    double operator()(const Tbd& b1d1) {
        bool isStrainNode = thrust::get<0>(b1d1);
        if (isStrainNode)
            return thrust::get<1>(b1d1);
        else
            return 0.0;
    }
};

// Functor for calculating the norm (magnitude) of a CVec3 tuple.
struct NormFunctor {
    __host__ __device__
    double operator()(const CVec3& vec) {
        // Calculate the Euclidean distance of the CVec3 tuple from the origin (0,0,0).
        // The magnitude of the CVec3 is returned as the norm.
        double result = sqrt(
            thrust::get<0>(vec) * thrust::get<0>(vec) +
            thrust::get<1>(vec) * thrust::get<1>(vec) +
            thrust::get<2>(vec) * thrust::get<2>(vec)
        );
        return result;
    }
};

//
//This part of the code includes functors for element-wise addition, vector norms, and random number generation. Each functor defines its operation and data type, enabling them to be used with Thrust algorithms like `thrust::transform` or `thrust::reduce`.
//
//The functors' functionalities are as follows:
//
//1. `UCVec3Add`: Takes two UCVec3 tuples and returns their element-wise sum as another UCVec3 tuple.
//2. `CVec3Add`: Takes two CVec3 tuples and returns their element-wise sum as another CVec3 tuple.
//3. `CVec4Add`: Takes two CVec4 tuples and returns their element-wise sum as another CVec4 tuple.
//4. `CVec3NormBinary`: Calculates the norm (magnitude) of a CVec3 tuple as a binary function. It computes the Euclidean distance between two CVec3 tuples and returns the magnitude.
//5. `CVec3NormUnary`: Calculates the norm (magnitude) of a CVec3 tuple as a unary function. It computes the Euclidean distance of the CVec3 tuple from the origin (0,0,0) and returns the magnitude.
//6. `psrnormgen`: Functor for generating normally distributed random numbers within a specified range.
//7. `psrunifgen`: Functor for generating uniformly distributed random numbers within a specified range.
//8. `TorsionAngleFunctor`: Calculates the torsion angle of a triplet of nodes. It takes the X, Y, and Z coordinates of three nodes and returns the torsion angle (radians) between them.
//9. `AveStrainFunctor`: Functor for averaging strain values. It takes a bool (indicating whether the node is a strain node) and a double (strain value) and returns the strain value if the node is a strain node; otherwise, it returns 0.0.
//10. `NormFunctor`: Calculates the norm (magnitude) of a CVec3 tuple. It computes the Euclidean distance of the CVec3 tuple from the origin (0,0,0) and returns the magnitude.
//
//These functors are essential building blocks used in various parts of the code to perform vector operations, random number generation, and angle calculations. They enable efficient parallel processing using Thrust's algorithms and can be applied to GPU-accelerated computations.


// Functor for checking if a value is greater than a given level.
struct IsGreaterThanLevel {
    double limit;

    __host__ __device__
    IsGreaterThanLevel(double& _limit) : limit(_limit) {}

    __device__
    bool operator()(double zPos) {
        // Returns true if zPos is greater than the specified limit.
        return (zPos > limit);
    }
};

// Functor for checking if a value is less than a given level.
struct IsLessThanLevel {
    double limit;

    __host__ __device__
    IsLessThanLevel(double& _limit) : limit(_limit) {}

    __device__
    bool operator()(double zPos) {
        // Returns true if zPos is less than the specified limit.
        return (zPos < limit);
    }
};

// Functor for checking if two Tuu tuples are equal.
struct tupleEqual {
    __host__ __device__
    bool operator()(Tuu x, Tuu y) {
        // Returns true if both elements of the two tuples are equal.
        return ((thrust::get<0>(x) == thrust::get<0>(y)) && (thrust::get<1>(x) == thrust::get<1>(y)));
    }
};

// Functor for checking if a value is equal to one.
struct IsEqualToOne {
    __host__ __device__
    bool operator()(const int& x) {
        // Returns true if x is not equal to one.
        return (x != 1);
    }
};

// Functor for checking if a value is not equal to zero.
struct isNotEqualZero {
    __host__ __device__
    bool operator()(const int& x) {
        // Returns true if x is not equal to zero.
        return (x != 0);
    }
};

// Functor for checking if a value is equal to zero.
struct isEqualZero {
    __host__ __device__
    bool operator()(const int& x) {
        // Returns true if x is equal to zero.
        return (x == 0);
    }
};

// Functor for checking if a value is greater than a given limit.
struct is_greater_than {
    int limit;

    __host__ __device__
    is_greater_than(int& _limit) : limit(_limit) {}

    __device__
    bool operator()(const int& x) {
        // Returns true if x is greater than the specified limit.
        if (x > limit) {
            return true;
        } else {
            return false;
        }
    }
};

// Functor for checking if a value is less than a given limit.
struct is_less_than {
    int limit;

    __host__ __device__
    is_less_than(int& _limit) : limit(_limit) {}

    __device__
    bool operator()(const int& x) {
        // Returns true if x is less than the specified limit.
        if (x < limit) {
            return true;
        } else {
            return false;
        }
    }
};

// Functor for calculating the inner product of a CVec3 tuple with itself.
struct CVec3InnerProduct {
    __host__ __device__
    double operator()(CVec3& v1) {
        // Returns the inner product of the CVec3 tuple with itself.
        return (thrust::get<0>(v1) * thrust::get<0>(v1) +
                thrust::get<1>(v1) * thrust::get<1>(v1) +
                thrust::get<2>(v1) * thrust::get<2>(v1));
    }
};

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
// FUNCTIONS USED IN AREA TRIANGLES AND BENDING TRIANGLES.

// Inline function for computing the cross product of two CVec3 tuples.
__host__ __device__
inline CVec3 CVec3_cross(CVec3 v1, CVec3 v2) {
    return thrust::make_tuple(
        thrust::get<1>(v1) * thrust::get<2>(v2) - thrust::get<2>(v1) * thrust::get<1>(v2),
        -(thrust::get<0>(v1) * thrust::get<2>(v2) - thrust::get<2>(v1) * thrust::get<0>(v2)),
        thrust::get<0>(v1) * thrust::get<1>(v2) - thrust::get<1>(v1) * thrust::get<0>(v2)
    );
};

// Inline function for element-wise multiplication of two CVec3 tuples.
__host__ __device__
inline CVec3 CVec3_mult(CVec3 v1, CVec3 v2) {
    return thrust::make_tuple(
        thrust::get<0>(v1) * thrust::get<0>(v2),
        thrust::get<1>(v1) * thrust::get<1>(v2),
        thrust::get<2>(v1) * thrust::get<2>(v2)
    );
};

// Inline function for scalar multiplication of a CVec3 tuple with a constant.
__host__ __device__
inline CVec3 CVec3_scalermult(double c, CVec3 v2) {
    return thrust::make_tuple(
        c * thrust::get<0>(v2),
        c * thrust::get<1>(v2),
        c * thrust::get<2>(v2)
    );
};

// Inline function for element-wise division of two CVec3 tuples.
__host__ __device__
inline CVec3 CVec3_div(CVec3 v1, CVec3 v2) {
    return thrust::make_tuple(
        thrust::get<0>(v1) / thrust::get<0>(v2),
        thrust::get<1>(v1) / thrust::get<1>(v2),
        thrust::get<2>(v1) / thrust::get<2>(v2)
    );
};

// Inline function for computing the dot product of two CVec3 tuples.
__host__ __device__
inline double CVec3_dot(CVec3 v1, CVec3 v2) {
    return (
        thrust::get<0>(v1) * thrust::get<0>(v2) +
        thrust::get<1>(v1) * thrust::get<1>(v2) +
        thrust::get<2>(v1) * thrust::get<2>(v2)
    );
};

// Inline function for element-wise addition of two CVec3 tuples.
__host__ __device__
inline CVec3 CVec3_plus(CVec3 v1, CVec3 v2) {
    return thrust::make_tuple(
        thrust::get<0>(v1) + thrust::get<0>(v2),
        thrust::get<1>(v1) + thrust::get<1>(v2),
        thrust::get<2>(v1) + thrust::get<2>(v2)
    );
}; 

// Inline function for element-wise addition of a constant to a CVec3 tuple.
__host__ __device__
inline CVec3 CVec3_plus(CVec3 v1, double d) {
    return thrust::make_tuple(
        thrust::get<0>(v1) + d,
        thrust::get<1>(v1) + d,
        thrust::get<2>(v1) + d
    );
};

// Inline function for element-wise addition of three CVec3 tuples.
__host__ __device__
inline CVec3 CVec3_plus(CVec3 v1, CVec3 v2, CVec3 v3) {
    return thrust::make_tuple(
        thrust::get<0>(v1) + thrust::get<0>(v2) + thrust::get<0>(v3),
        thrust::get<1>(v1) + thrust::get<1>(v2) + thrust::get<1>(v3),
        thrust::get<2>(v1) + thrust::get<2>(v2) + thrust::get<2>(v3)
    );
};

// Inline function for element-wise addition of six CVec3 tuples.
__host__ __device__
inline CVec3 CVec3_plus(CVec3 v1, CVec3 v2, CVec3 v3, CVec3 v4, CVec3 v5, CVec3 v6) {
    return thrust::make_tuple(
        thrust::get<0>(v1) + thrust::get<0>(v2) + thrust::get<0>(v3) + thrust::get<0>(v4) + thrust::get<0>(v5) + thrust::get<0>(v6),
        thrust::get<1>(v1) + thrust::get<1>(v2) + thrust::get<1>(v3) + thrust::get<1>(v4) + thrust::get<1>(v5) + thrust::get<1>(v6),
        thrust::get<2>(v1) + thrust::get<2>(v2) + thrust::get<2>(v3) + thrust::get<2>(v4) + thrust::get<2>(v5) + thrust::get<2>(v6)
    );
};

// Inline function for element-wise subtraction of two CVec3 tuples.
__host__ __device__
inline CVec3 CVec3_subtract(CVec3 v1, CVec3 v2) {
    return thrust::make_tuple(
        thrust::get<0>(v1) - thrust::get<0>(v2),
        thrust::get<1>(v1) - thrust::get<1>(v2),
        thrust::get<2>(v1) - thrust::get<2>(v2)
    );
};

// Inline function for element-wise subtraction of two CVec3 tuples (alternative implementation).
__host__ __device__
inline CVec3 CVec3_minus(CVec3 v1, CVec3 v2) {
    return thrust::make_tuple(
        thrust::get<0>(v1) - thrust::get<0>(v2),
        thrust::get<1>(v1) - thrust::get<1>(v2),
        thrust::get<2>(v1) - thrust::get<2>(v2)
    );
};

#endif /* SYSTEMSTRUCTURES_H_ */
