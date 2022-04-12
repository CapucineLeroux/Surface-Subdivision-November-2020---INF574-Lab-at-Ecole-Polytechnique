/* ========================================================================= *
 *                                                                           *
 *                       Luca Castelli Aleardi                       		 *
 *           Copyright (c) 2019, Ecole Polytechnique                		 *
 *           Department of Computer Science                  				 *
 *                          All rights reserved.                             *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * This file is part of the course material developed for		             *
 *   INF574 Digital Representation and Analysis of Shapes (2019/20)			 *
 * ========================================================================= */
#include <igl/opengl/glfw/Viewer.h>

#ifndef HALFEDGE_DS_HEADER
  #define HALFEDGE_DS_HEADER
  #include "HalfedgeDS.cpp"
#endif

using namespace Eigen;
using namespace std;

/**
 * @author Luca Castelli Aleardi (2019)
 */
class LoopSubdivision
{

public:
    /** 
	 * Initialize the data structures
	 **/
    LoopSubdivision(MatrixXd &V_original, MatrixXi &F_original, HalfedgeDS &mesh)
    {
        he = &mesh;
        V = &V_original;
		F = &F_original; // NOT NEEDED if using the half-edge data structure
        int e = he->sizeOfHalfedges() / 2; // number of edges in the original mesh
        int n = V_original.rows();         // number of vertices in the original mesh

        V1 = new MatrixXd();
        V1->setZero(n+e,3);

        F1 = new MatrixXi();
        F1->setZero(4.*(he->sizeOfFaces()),3);

        nVertices = he->sizeOfVertices();
        nFaces = he->sizeOfFaces();

    }

    /** 
	 * Perform the subdivision of the mesh (just perform one round subdivision). <b>
	 * As result, the resulting subdivided mesh is stored in arrays 'V1' and 'F1' 
	 **/
    void subdivide()
    {
        std::cout << "Performing one round subdivision" << endl;
        int e = he->sizeOfHalfedges() / 2; // number of edges in the original mesh
        int n = he->sizeOfVertices();      // number of vertices in the original mesh
        int nF = he->sizeOfFaces();         // number of vertices in the original mesh

        // first step: perform a copy of the original points
        for (int i=0 ; i<n ; i++){
            V1->row(i) = V->row(i);
        }

        // second step: compute new midpoint vertices and assign a number, between 0..e-1, to all halfedges

        MatrixXi HE = MatrixXi::Zero(2*e,1); //will store the new midpoints indices in the corresponding pair of halfedges
        MatrixXd midpoint; //will be the new computed midpoint of the edge
        int h_opp;
        int cpt = 1; //count the number of edges updated (starts at one to avoid the edge case of the 1st halfedge)
        //edge case
        midpoint = computeEdgePoint(0);
        V1->row(n) = midpoint;
        //for each halfedge
        for (int h=0 ; h<2*e ; h++){

            h_opp = he->getOpposite(h);

            //if the pair of halfedges haven't been updated with a pointer to the new midpoint, we do so in HE
            //and we put the correponding midpoint at the right place in V1
            //takes the edge case of the 1st edge into account
            if (HE(h,0)==0 && HE(h_opp,0)==0 && h!=0 && h_opp!=0){

                //update pointers of the pair of halfedges
                HE(h,0) = cpt;
                HE(h_opp,0) = cpt;

                //put the midpoint in V1
                midpoint = computeEdgePoint(h);
                V1->row(n+cpt) = midpoint;

                cpt += 1;
            }
        }

        //udpate the original points
        for (int i=0 ; i<n ; i++){
            V1->row(i) = updateOriginalPoint(i);
        }


        // third step: set the face/vertex incidence relations
        int v0;
        int v1;
        int v2;
        int v3;
        int v4;
        int v5;
        int h0;
        int h1;
        int h2;
        for (int i=0 ; i<nF ; i++){
            //get the 3 vertices of the original ith face
            v0 = F->coeffRef(i,0);
            v1 = F->coeffRef(i,1);
            v2 = F->coeffRef(i,2);

            //get the 3 halfedges interior to the original ith face
            h0 = 3*i;
            h1 = 3*i+1;
            h2 = 3*i+2;

            //get the 3 new midpoints in the ith face
            v3 = n + HE(h0,0);
            v4 = n + HE(h1,0);
            v5 = n + HE(h2,0);

            //create the 4 new faces
            F1->row(4*i) << v1,v4,v3;
            F1->row(4*i+1) << v4,v2,v5;
            F1->row(4*i+2) << v3,v5,v0;
            F1->row(4*i+3) << v3,v4,v5;
        }

    }

    /** 
	 * Return the number of half-edges
	 **/
    MatrixXd getVertexCoordinates()
    {
        return *V1;
    }

    /** 
	 * Return the number of faces
	 **/
    MatrixXi getFaces()
    {
        return *F1;
    }

    /** 
	 * Print the combinatorial information of the subdivided mesh <b>
	 * verbosity=0: print only the number of vertices and faces <b>
	 * verbosity=1: print all incidence relations
	 **/
    void print(int verbosity)
    {
        cout << "\tn=" << nVertices << ", f=" << nFaces << endl;

        if (verbosity > 0) // print all vertex coordinates and face/vertex incidence relations
        {
            for (int i = 0; i < nVertices; i++)
            {
                cout << "v" << i << ": " << V1->row(i) << endl;
            }

            std::cout << "new faces: " << nFaces << endl;
            for (int i = 0; i < nFaces; i++)
            {
                cout << "f" << i << ": " << F1->row(i) << endl;
            }
        }
    }

private:
    /**
	 * Compute the midpoint of the given half-edge 'h=(u,v)'
	 */
    MatrixXd computeEdgePoint(int h)
    {
        //extract the 4 vertices needed to subdivide
        int v1 = he->getTarget(h);
        MatrixXd s1 = V->row(v1);

        int h2 = he->getOpposite(h);
        int v2 = he->getTarget(h2);
        MatrixXd s2 = V->row(v2);

        int h3 = he->getNext(h);
        int v3 = he->getTarget(h3);
        MatrixXd s3 = V->row(v3);

        int h4 = he->getNext(h2);
        int v4 = he->getTarget(h4);
        MatrixXd s4 = V->row(v4);

        //create the new midpoint of the edge
        MatrixXd s5 = (3./8.)*(s1+s2) + (1./8.)*(s3+s4);
        //float N = s5.norm();
        //s5 = s5/N;

        return s5;

    }

    /**
	 * Given a vertex 'v' of the original mesh, compute and return its new coordinates
	 */
    MatrixXd updateOriginalPoint(int v)
    {
        //get information about the former point indexed by v
        MatrixXd p = V->row(v); //former point coordinates
        int d = vertexDegree(v);
        int h = he->getEdge(v);
        int h_bis = he->getOpposite(h);

        //calculate alpha according to the degree
        float alpha;
        if (d==3){
            alpha = 3./16.;
        }
        else if (d>3){
            alpha = 3./(8.*d);
        }

        //find the new point p1
        int v1 = he->getTarget(h_bis);
        MatrixXd p1 = V->row(v1);//new point coordinates
        for (int i=1 ; i<d ; i++){
            h_bis = he->getNext(he->getOpposite(h_bis));
            p1 += V->row(he->getTarget(h_bis));
        }
        p1 = (1.-alpha*d)*p + alpha*p1;

        return p1;
    }

    int vertexDegree(int v)
    {
        int result = 0;

        int h = he->getEdge(v);
        int h_bis = he->getOpposite(he->getNext(h));

        while (h_bis != h){
            h_bis = he->getOpposite(he->getNext(h_bis));
            result++;
        }

        return result+1;

    }

    /** Half-edge representation of the original input mesh */
    HalfedgeDS *he;
    MatrixXd *V; // vertex coordinates of the original input mesh

	/** faces/vertex incidence relations in the original mesh */
	MatrixXi *F; // REMARK: not needed if using the half-edge data structure

    int nVertices, nFaces; // number of vertices, faces in the new subdivided mesh
    MatrixXd *V1;          // vertex coordinates of the new subdivided mesh
    MatrixXi *F1;          // faces of the new subdivided mesh
};
