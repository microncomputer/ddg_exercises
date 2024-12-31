// Implement member functions for SimplicialComplexOperators class.
#include "simplicial-complex-operators.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

/*
 * Assign a unique index to each vertex, edge, and face of a mesh.
 * All elements are 0-indexed.
 *
 * Input: None. Access geometry via the member variable <geometry>, and pointer to the mesh via <mesh>.
 * Returns: None.
 */
void SimplicialComplexOperators::assignElementIndices() {

    // Needed to access geometry->vertexIndices, etc. as cached quantities.
    // Not needed if you're just using v->getIndex(), etc.
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    // You can set the index field of a vertex via geometry->vertexIndices[v], where v is a Vertex object (or an
    // integer). Similarly you can do edges and faces via geometry->edgeIndices, geometry->faceIndices, like so:
    size_t idx = 0;
    for (Vertex v : mesh->vertices()) {
        idx = geometry->vertexIndices[v];
    }

    for (Edge e : mesh->edges()) {
        idx = geometry->edgeIndices[e];
    }

    for (Face f : mesh->faces()) {
        idx = geometry->faceIndices[f];
    }

    // You can more easily get the indices of mesh elements using the function getIndex(), albeit less efficiently and
    // technically less safe (although you don't need to worry about it), like so:
    //
    //      v.getIndex()
    //
    // where v can be a Vertex, Edge, Face, Halfedge, etc. For example:

    for (Vertex v : mesh->vertices()) {
        idx = v.getIndex(); // == geometry->vertexIndices[v])
    }

    // Geometry Central already sets the indices for us, though, so this function is just here for demonstration.
    // You don't have to do anything :)
}

/*
 * Construct the unsigned vertex-edge adjacency matrix A0.
 *
 * Input:
 * Returns: The sparse vertex-edge adjacency matrix which gets stored in the global variable A0.
 * note: the naming is inconsistent with the function that follows this one. I believe it would be more correct to call this the edge-vertex matrix 
 * since it produces the |E|x|V| adjacency matrix.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildVertexEdgeAdjacencyMatrix() const {

    // TODO
    // Note: You can build an Eigen sparse matrix from triplets, then return it as a Geometry Central SparseMatrix.
    // See <https://eigen.tuxfamily.org/dox/group__TutorialSparse.html> for documentation.

    // renaming the triplet type to T for simplicity
    typedef Eigen::Triplet<size_t> T; 
    std::vector<T> tripletList;
    
    // get num vertices and num edges and make sparse mat of that size
    SparseMatrix<size_t> spMat(mesh->nEdges(), mesh->nVertices());
    
    // for each edge look at adjacent vertices and add to that row
    for(Edge e: mesh->edges()){
        for(Vertex v: e.adjacentVertices()){
            tripletList.push_back(T(e.getIndex(), v.getIndex(), 1));
        }
    }

    spMat.setFromTriplets(tripletList.begin(), tripletList.end());

    return spMat;
}

/*
 * Construct the unsigned face-edge adjacency matrix A1.
 *
 * Input:
 * Returns: The sparse face-edge adjacency matrix which gets stored in the global variable A1.
 * 
 * Note: The name of this function is misleading with respect to the last one. buildVertexEdge makes an ExV matrix but buildFaceEdge does not build an ExF matrix. It builds an FxE. 
 * in that way.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildFaceEdgeAdjacencyMatrix() const {
    typedef Eigen::Triplet<size_t> T; 
    std::vector<T> tripletList;
    //tripletList.reserve(mesh->nFaces()*mesh->nEdges());
    
    // get num vertices and num edges and make sparse mat of that size
    SparseMatrix<size_t> spMat(mesh->nFaces(), mesh->nEdges());
    
    // for each face look at adjacent edges and add to that row
    for(Face f: mesh->faces()){
        for(Edge e: f.adjacentEdges()){
            tripletList.push_back(T(f.getIndex(), e.getIndex(), 1));
        }
    }

    spMat.setFromTriplets(tripletList.begin(), tripletList.end());
    
    return spMat;
}

/*
 * Construct a vector encoding the vertices in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |V|, where |V| = # of vertices in the mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildVertexVector(const MeshSubset& subset) const {

    Vector<size_t> vertices = Vector<size_t>::Zero(mesh->nVertices());
    for(int idx: subset.vertices){
        vertices(idx) = 1;
    }
    return vertices;
}

/*
 * Construct a vector encoding the edges in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |E|, where |E| = # of edges in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildEdgeVector(const MeshSubset& subset) const {

    Vector<size_t> edges = Vector<size_t>::Zero(mesh->nEdges());
    for(int idx: subset.edges){
        edges(idx) = 1;
    }
    return edges;
}

/*
 * Construct a vector encoding the faces in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |F|, where |F| = # of faces in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildFaceVector(const MeshSubset& subset) const {
    
    Vector<size_t> faces = Vector<size_t>::Zero(mesh->nFaces());
    for(int idx: subset.faces){
        faces(idx) = 1;
    }
    return faces;
}

/*
 * Compute the simplicial star St(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The star of the given subset.
 */
MeshSubset SimplicialComplexOperators::star(const MeshSubset& subset) const {
    // the star(subset) is the collection of all simplices of the complex (ie. the mesh) that contain ANY simplex in subset
   
    MeshSubset Star;

    // add each vertex in subset
    Vector<size_t> vertices = buildVertexVector(subset);
    for(size_t i = 0; i <mesh->nVertices(); ++i){
        if (vertices[i]) Star.addVertex(i);
    }

    // for each edge in subset, add them and adjacent faces.
    Vector<size_t> edges = buildEdgeVector(subset) + A0 * vertices;
    for(size_t i = 0; i <mesh->nEdges(); ++i){
        if (edges[i]) Star.addEdge(i);
    }

    // add the rest of the faces in subset.
    Vector<size_t> faces = buildFaceVector(subset) + A1 * edges;
    for(size_t i = 0; i <mesh->nFaces(); ++i){
        if (faces[i]) Star.addFace(i);
    }
    
    return Star; 
}


/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the given subset.
 */
MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {
    //start by adding everything in subset to the closure
    MeshSubset Closure = subset.deepCopy();
    
    Vector<size_t> vertices = buildVertexVector(subset);
    Vector<size_t> edges = buildEdgeVector(subset);
    Vector<size_t> faces = buildFaceVector(subset);


    // for each face, add the adjacent edges
    edges = edges + (faces.transpose() * A1).transpose();
    for(size_t i= 0; i < edges.size(); ++i){
        if(edges[i])
            Closure.addEdge(i);
    }

    // for each edge, add the adjacent vertices
    vertices = vertices + (edges.transpose() * A0).transpose();
    for(size_t i= 0; i < vertices.size(); ++i){
        if(vertices[i])
            Closure.addVertex(i);
    }

    return Closure; // placeholder
}

/*
 * Compute the link Lk(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The link of the given subset.
 */
MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {
    //link = closure(star(subset)) setminus star(closure(subset))
    MeshSubset clofst = closure(star(subset));
    MeshSubset stofcl = star(closure(subset));
    clofst.deleteSubset(stofcl);
    return clofst;
}

/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: True if given subset is a simplicial complex, false otherwise.
 */
bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {
    if(subset.equals(closure(subset)))
        return true;
    return false; 
}

/*
 * Check if the given subset S is a pure simplicial complex. If so, return the degree of the complex. Otherwise, return
 * -1.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: int representing the degree of the given complex (-1 if not pure)
 */
int SimplicialComplexOperators::isPureComplex(const MeshSubset& subset) const {
    //check that it is a complex
    if(!isComplex(subset))
        return -1;
    
    Vector<size_t> vertices = buildVertexVector(subset);
    Vector<size_t> edges = buildEdgeVector(subset);
    Vector<size_t> faces = buildFaceVector(subset);

    int degree=0;
    //if there are any faces, it will need to be a pure 2-complex
    if(subset.faces.size() > 0){
         degree = 2;
        //check that each edge in subset is adjacent to a face in subset
        for(size_t e: subset.edges){
            //if col e of A1 is empty, there's no face connecting edge e at all
            if(!A1.innerVector(e).nonZeros()) return -1;
        }
        Vector<size_t> edgesadjfacesinsubset = (faces.transpose() * A1).transpose();
        //if edges adjacent to faces in subset is not equal to subset's EdgeVector, it is not pure
        if(!(edgesadjfacesinsubset.count() == subset.edges.size())) return -1;
    }
    //else if there are no faces but any edges, it will need to be a pure 1-complex
    else if(subset.edges.size() > 0)
        degree = 1;
    else degree = 0;
    if(degree == 1 || degree == 2){
        //check that each vertex in subset is adjacent to an edge in subset
        for(size_t v: subset.vertices){
            //if col v of A0 is empty, there's no edge connecting vertex v at all
            if(!A0.innerVector(v).nonZeros()) return -1;
        }
        Vector<size_t> verticesadjedgesinsubset = (edges.transpose() * A0).transpose();
        //if vertices adjacent to edges in subset is not equal to subset's VertexVector, it is not pure
        if(!(verticesadjedgesinsubset.count() == subset.vertices.size())) return -1;
    }
    
    //else it is pure 0-complex
    return degree; 
}

/*
 * Compute the set of simplices contained in the boundary bd(S) of the selected subset S of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The boundary of the given subset.
 */
MeshSubset SimplicialComplexOperators::boundary(const MeshSubset& subset) const {
    MeshSubset boundary;
    int k = isPureComplex(subset);
    if (k == -1 || k == 0) return boundary;
    
    Vector<size_t> edges = buildEdgeVector(subset);

    //depending on degree k, look at the respective adj mat to get the (k-1)-simplices in the boundary.
    if(k == 1){
        //for each vertex in vertices, if there is more than one adjacent edge, it is not in the boundary. else it is.
        Vector<size_t> vertices = buildVertexVector(subset);
        SparseMatrix<size_t> subsetadjmat = (edges.asDiagonal())*A0;
        subsetadjmat = subsetadjmat *(vertices.asDiagonal());
        subsetadjmat = subsetadjmat.pruned();
        for(size_t v: subset.vertices){
            //if there is just one adjacent edge, it is a boundary vertex 
            if(subsetadjmat.innerVector(v).nonZeros()==1)
                boundary.addVertex(v);
        }

    }
    else if(k == 2){
        //for each edge in edges, if there is more than one adjacent face, it is not in the boundary. else it is.
        Vector<size_t> faces = buildFaceVector(subset);
        SparseMatrix<size_t> subsetadjmat = faces.asDiagonal()*A1*edges.asDiagonal();
        subsetadjmat = subsetadjmat.pruned();
        for(size_t e: subset.edges){
            //if there is just one adjacent face, it is a boundary edge 
            if(subsetadjmat.innerVector(e).nonZeros()==1){
                boundary.addEdge(e);
                //and add the vertices of each boundary edge
                for (Vertex v: mesh->edge(e).adjacentVertices()){
                    boundary.addVertex(geometry->vertexIndices[v]);
                }
            }
        }
    }
    return boundary; // placeholder
}