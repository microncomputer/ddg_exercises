// PLEASE READ:
//
// This file additional geometry routines for the VertexPositionGeometry class in Geometry Central. Because we are
// "inside" the class, we no longer have to call
//
//          geometry->inputVertexPositions[v], etc.
//
// We can just call
//
//          this->inputVertexPositions[v], etc.
//
// or simply
//
//          inputVertexPositions[v], etc.
//
// In addition, we no longer access the corresponding surface mesh via
//
//          mesh->vertices(), etc.
//
// but instead <mesh> is not a pointer anymore, so we use
//
//          mesh.vertices(), etc.
//
// Functions in this file can be called from other projects simply by using geometry->buildHodgeStar0Form(), etc. where
// "geometry" is a pointer to a VertexPositionGeometry. This avoids having to declare a GeometryRoutines object in every
// project, and also mimics the way that geometry routines are normally called in Geometry Central.
//
// Other notes: In this file, you can use the constant pi by using PI.

#include "geometrycentral/surface/vertex_position_geometry.h"


namespace geometrycentral {
namespace surface {


/*
 * Build Hodge operator on 0-forms.
 * By convention, the area of a vertex is 1.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar0Form() const {
    // compute regions of the dual faces to each vertex and 
    //SparseMatrix<double> H0(mesh.nVertices(),1);
    SparseMatrix<double> H0(mesh.nVertices(),mesh.nVertices());
    std::vector<Eigen::Triplet<double>> tripletlist;
    for (Vertex v: mesh.vertices()){
        tripletlist.push_back(Eigen::Triplet<double>(v.getIndex(),v.getIndex(), barycentricDualArea(v)));
    }
    H0.setFromTriplets(tripletlist.begin(), tripletlist.end());
    return H0;
}

/*
 * Build Hodge operator on 1-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar1Form() const {
    SparseMatrix<double> H1(mesh.nEdges(),mesh.nEdges());

    std::vector<Eigen::Triplet<double>> tripletlist;
    for (Edge e: mesh.edges()){
        //found edgecotanweight in the source code for vertexpositiongeometry 
        double cotanweight_e = edgeCotanWeight(e);
        tripletlist.push_back(Eigen::Triplet<double>(e.getIndex(),e.getIndex(), cotanweight_e));
    }
    H1.setFromTriplets(tripletlist.begin(), tripletlist.end());
    return H1;   
}

/*
 * Build Hodge operator on 2-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 2-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar2Form() const {
    SparseMatrix<double> H2(mesh.nFaces(),mesh.nFaces());

    std::vector<Eigen::Triplet<double>> tripletlist;
    for (Face f: mesh.faces()){
        // size of dual face (a vertex) over size of face
        // 1/2 (sum|ei|)
        double face_area = faceArea(f);
        /*
        double dualarea = 0;
        for(Edge e: f.adjacentEdges()){
            dualarea += edgeLength(e);
        }
        dualarea = dualarea/2;
        */

        tripletlist.push_back(Eigen::Triplet<double>(f.getIndex(),f.getIndex(), 1/face_area));
    }
    H2.setFromTriplets(tripletlist.begin(), tripletlist.end());
    
    return H2;  
}

/*
 * Build exterior derivative on 0-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative0Form() const {
 
    // for a discrete 0-form phi (values at vertices), the derivative is a 1-form. 
    // it is the integrals of dphi along edges.
    // ie. add up the values of applying the input differential form to the k-form
    // returns a k+1 form, 1-form in this case
    /* steps:
    1. make a ExV matrix 
    */
    
    //it's ExV because it is mapping values from V to values from E 
    SparseMatrix<double> d0(mesh.nEdges(),mesh.nVertices());
    std::vector<Eigen::Triplet<double>> tripletlist;

    for(Edge e: mesh.edges()){
        tripletlist.push_back(Eigen::Triplet<double>(e.getIndex(), e.halfedge().tipVertex().getIndex(), 1));
        tripletlist.push_back(Eigen::Triplet<double>(e.getIndex(), e.halfedge().tailVertex().getIndex(), -1));
    }
    d0.setFromTriplets(tripletlist.begin(), tripletlist.end());
    return d0;
}

/*
 * Build exterior derivative on 1-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative1Form() const {
    SparseMatrix<double> d1(mesh.nFaces(),mesh.nEdges());
    std::vector<Eigen::Triplet<double>> tripletlist;

    for(Face f: mesh.faces()){
        for(Edge e: f.adjacentEdges()){
            if(e.halfedge().face() == f) //it's stored as 1 or 0
                tripletlist.push_back(Eigen::Triplet<double>(f.getIndex(),e.getIndex(), 1));
            else
                tripletlist.push_back(Eigen::Triplet<double>(f.getIndex(),e.getIndex(), -1));
        }
    }
    d1.setFromTriplets(tripletlist.begin(), tripletlist.end());
    return d1;
}

} // namespace surface
} // namespace geometrycentral