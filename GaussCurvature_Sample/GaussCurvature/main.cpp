/*!
*      \file main.cpp
*      \brief Main program for Slit Mapping/Riemann Mapping/Harmonic Mapping
*	   \author David Gu
*      \date Document 10/07/2010
*
*		The code performs the following tasks
*
*       1. Compute the harmonic maps for simply connected domain with a single boundary component to the unit disk.
*/

/*!	Harmonic Mapping
 */
#include "Conformal/GaussCurvature/GaussCurvatureMesh.h" //Harmonic Mapping
#include "Conformal/GaussCurvature/GaussCurvature.h" //Harmonic Mapping

#include "Mesh/iterators.h"
#include <vector>
#include <list>
#include <stdlib.h>
#include <math.h>
using namespace MeshLib;

typedef CBaseMesh<CVertex, CEdge, CFace, CHalfEdge> delaunay;
//typedef FaceVertexIterator<M> FaceVertexIterator;
typedef CGaussCurvatureMesh<CVertex,CEdge,CFace,CHalfEdge> M;

delaunay Delaunay;

int id_f = 0;
int id_v = 0;

bool LegalizeEdge(CVertex*, CEdge*);
double computeS(CPoint2, CPoint2, CPoint2);

/*!	\brief Helper function to remind the users of the usage of the commands	
 *
 */

void help(char * exe )
{
	printf("Usage:\n");
	printf("%s -gauss_curvature  input_mesh \n", exe );
}

//compute the total Gauss curvature
void _Gauss_Curvature( const char * _input)
{
	CGCMesh mesh;
	mesh.read_m( _input );

	CGaussCurvature<CGCMesh> mapper( & mesh );

	mapper._calculate_face_normal();
	mapper._calculate_vertex_normal();
	mapper._calculate_Euler_characteristics();
	mapper._calculate_curvature();
	//system("pause");
};

CFace* createFace1(std::vector<CVertex*> v)
{
	CFace* f;
	f = Delaunay.createFace(v,id_f);
	id_f++;
	return f;
}

CVertex* createVertex1(CPoint2 pv)
{
	CVertex* v;
    v = Delaunay.createVertex(id_v);
	v->uv()=pv;
	id_v++;
	return v;
}

CFace* LocatePoint(CVertex* v){
	CPoint2 p = v->uv();
	CFace* f;
	std::list<CFace*> faces;
	faces.push_back(Delaunay.faces().back());
	//std::cout<<faces.back()<<std::endl;
	std::vector<CVertex*> vector1;
	do{
		f = faces.back();
		faces.pop_back();
		for( FaceVertexIterator<M> viter(f); !viter.end(); ++ viter){
			vector1.push_back(*viter);
		}
		CVertex *v0,*v1,*v2;
		CPoint2 p0,p1,p2;
		v0 = vector1[0];
		v1 = vector1[1];
		v2 = vector1[2];
		p0 = v0->uv();
		p1 = v1->uv();
		p2 = v2->uv();
		vector1.clear();

		double result0,result1,result2,uparea;
		double area = computeS(p0,p1,p2);

		uparea = computeS(p,p0,p1);
		result0 = uparea/area;
		if (result0 < 0){
			CFace *newface;

			CHalfEdge *he;

			he = Delaunay.vertexHalfedge(v0,v1);
			if(he->he_next()->target() == v2)
				newface = he->he_sym()->face();
			else
				newface = he->face();
			faces.push_back(newface);
		}
		else{
			uparea = computeS(p,p1,p2);
			result1 = uparea/area;
			if (result1 < 0){
				CFace *newface;

				CHalfEdge *he;

				he = Delaunay.vertexHalfedge(v1,v2);
				if(he->he_next()->target() == v0)
					newface = he->he_sym()->face();
				else
					newface = he->face();
				faces.push_back(newface);

			}
			else{
				uparea = computeS(p,p2,p0);
				result2 = uparea/area;
				if (result2 < 0){
					CFace *newface;

					CHalfEdge *he;

					he = Delaunay.vertexHalfedge(v2,v0);
					if(he->he_next()->target() == v1)
						newface = he->he_sym()->face();
					else
						newface = he->face();
					faces.push_back(newface);
				}

				else if (result0 > 0 && result1 > 0 && result2 > 0)
					return f;
			}
		}

	}while(!faces.empty());

	return NULL;
}

void FaceSplit(CFace* pFace, CVertex* pv)
{
	std::vector<CVertex* > vec1;
	std::vector<CVertex* > vec2;
	std::vector<CVertex* > vec3;
	std::vector<CVertex* > v;
	for(FaceVertexIterator<M> fviter(pFace);!fviter.end();++ fviter)
	{
		CVertex *pV = *fviter;
		//CPoint2 p = pV->uv();
		//v.push_back(p);
		v.push_back(pV);
	}
	CVertex *pv0=v[0]; CVertex *pv1=v[1]; CVertex *pv2=v[2];	
	v.clear();
	//first delete the current face
	Delaunay.deleteFace(pFace);
	//then split the face into three faces
	vec1.clear();
	vec1.push_back(pv0);
	vec1.push_back(pv1);
	vec1.push_back(pv);
	CFace *face1 = createFace1(vec1);

	vec2.clear();
	vec2.push_back(pv1);
	vec2.push_back(pv2);
	vec2.push_back(pv);
	CFace *face2 = createFace1(vec2);

	vec3.clear();
	vec3.push_back(pv0);
	vec3.push_back(pv2);
	vec3.push_back(pv);
    CFace *face3 = createFace1(vec3);
    
	CEdge *e1 = Delaunay.vertexEdge(pv0,pv1);
	LegalizeEdge(pv,e1);
		
	CEdge *e2 = Delaunay.vertexEdge(pv1,pv2);
	LegalizeEdge(pv,e2);
	
	CEdge *e3 = Delaunay.vertexEdge(pv0,pv2);
	LegalizeEdge(pv,e3);
}

CEdge* edgeSwap(CEdge* e, CVertex* pv, CVertex* v2)
{
	CEdge* newEdge = NULL;
	if(Delaunay.isBoundary(e))
		newEdge = e;
	else{
	CFace *currentFace1;
	CFace *currentFace2;
	CFace *newFace1;
	CFace *newFace2;
	std::vector<CVertex*> vec1;
    std::vector<CVertex*> vec2;
	CVertex *v0 = Delaunay.edgeVertex1(e);
	CVertex *v1 = Delaunay.edgeVertex2(e);

	//find Vertex v2 according to the input edge
	
	
	CPoint2 p0 = v0->uv();
	CPoint2 p1 = v1->uv();
	CPoint2 p2 = v2->uv();

	//first delete the current two faces
	

	CHalfEdge *he = Delaunay.edgeHalfedge(e,0);
	CHalfEdge *he1 = Delaunay.edgeHalfedge(e,1);

//	currentFace1 = he->face();
//	currentFace2 = he1->face();
	currentFace1 = Delaunay.halfedgeFace(he);
	currentFace2 = Delaunay.halfedgeFace(he1);

	Delaunay.deleteFace(currentFace1);
	Delaunay.deleteFace(currentFace2);
	
	// then create two new faces
	vec1.push_back(pv);
    vec1.push_back(v0);
    vec1.push_back(v2);
	newFace1 = createFace1(vec1);

	vec2.push_back(pv);
    vec2.push_back(v1);
    vec2.push_back(v2);
	newFace2 = createFace1(vec2);

	newEdge = Delaunay.vertexEdge(pv,v2);
	}
	return newEdge;
}

double computeS(CPoint2 a,CPoint2 b,CPoint2 c)
{
	//CPoint2 a,b,c;
	double S;
	S = (a[0]-c[0])*(b[1]-c[1])-(a[1]-c[1])*(b[0]-c[0]);
	return S;
}

bool LegalizeEdge(CVertex* pv, CEdge* e)
{  
	// first decide whether the edge is on the boundary
	if(Delaunay.isBoundary(e))
		return false;
	else{
	CVertex *v0 = Delaunay.edgeVertex1(e);
	CVertex *v1 = Delaunay.edgeVertex2(e);
	CVertex* v2;
	//find Vertex v2 according to the input edge
	CHalfEdge *he = Delaunay.vertexHalfedge(v0, v1);
	CHalfEdge *nexthe = he->he_next();
	if(nexthe->target()==pv)
	v2 = he->he_sym()->he_next()->target();
	else
	v2 = nexthe->target();

	CEdge *e1 = Delaunay.vertexEdge(v1,v2);
	CEdge *e2 = Delaunay.vertexEdge(v0,v2);

	CPoint2 p0 = v0->uv();
	CPoint2 p1 = v1->uv();
	CPoint2 p2 = v2->uv();
	CPoint2 ppv = pv->uv();
/*
	CPoint2 ps0, ps1, ps2, ps3, ps4, ps5;*/

	CPoint2 ps0(ppv[0]*ppv[0]+ppv[1]*ppv[1],ppv[1]);
	CPoint2 ps1(p0[0]*p0[0]+p0[1]*p0[1],p0[1]);
	CPoint2 ps2(p1[0]*p1[0]+p1[1]*p1[1],p1[1]);
	CPoint2 ps3(ppv[0],ppv[0]*ppv[0]+ppv[1]*ppv[1]);
	CPoint2 ps4(p0[0],p0[0]*p0[0]+p0[1]*p0[1]);
	CPoint2 ps5(p1[0],p1[0]*p1[0]+p1[1]*p1[1]);

	//double S;
    double circlex = 0.0, circley = 0.0;
	double radius,distance;
	CPoint2 circle(circlex,circley);
	double S1, S2, S3;
    S1= computeS(ps0,ps1,ps2);
	S2= computeS(ppv,p0,p1);
	S3= computeS(ps3,ps4,ps5);
	circlex = S1/(2*S2);
	circley = 1.0/2*(S3/S2);
	radius = sqrt((ppv[0]-circlex)*(ppv[0]-circlex)+(ppv[1]-circley)*(ppv[1]-circley));
	distance = sqrt((p2[0]-circlex)*(p2[0]-circlex)+(p2[1]-circley)*(p2[1]-circley)); 
	if(distance>radius)
		return false;
	else{
		CEdge* newEdge = edgeSwap(e, pv, v2);
		LegalizeEdge(pv,e1);
		LegalizeEdge(pv,e2);
		return true;
	}
	}	
}

void InsertVertex(CVertex* pv)
{
	CFace* pFace = LocatePoint(pv);
	FaceSplit(pFace, pv);
	//Legalize three edges of the original face in FaceSplit Method.

}

void DelaunayTriangulation()
{
 //create a large triangle
// create 3 triangle faces
 double x0=0.0, y0=0.0, x1=50.0, y1=100.0, x2=100.0, y2=0.0, x3=50.0, y3=50.0;
 double a,b;

 CPoint2 p0(x0,y0);
 CPoint2 p1(x1,y1);
 CPoint2 p2(x2,y2);
 CPoint2 p3(x3,y3);
 
 std::vector<CVertex*> vec1;
 std::vector<CVertex*> vec2;
 std::vector<CVertex*> vec3;

 CVertex *v0 = createVertex1(p0);
 //v0->uv()= p0;

 CVertex *v1 = createVertex1(p1);
 //v1->uv()= p1;
 
 CVertex *v2 = createVertex1(p2);
 //v2->uv()= p2;

 CVertex *v3 = createVertex1(p3);
 //v3->uv() = p3;

 vec1.push_back(v0);
 vec1.push_back(v1);
 vec1.push_back(v3);
 CFace *face1 = createFace1(vec1);

 vec2.push_back(v1);
 vec2.push_back(v2);
 vec2.push_back(v3);
 CFace *face2 = createFace1(vec2);

 vec3.push_back(v0);
 vec3.push_back(v2);
 vec3.push_back(v3);
 CFace *face3 = createFace1(vec3);

 // repeat step 2 and 3
 int num = 5;
 for(int i=0;i<num;i++){
 //Randomly generate a point p in unit square
 a = 50+rand()/double(RAND_MAX);
 b = rand()/double(RAND_MAX);
  //a = 50.5;
  //b = 0.5;
 CPoint2 p_random(a,b);
 CVertex* pv = createVertex1(p_random);

// pv->uv() = p_random;
 
 //InsertVertex P
 InsertVertex(pv);
 }
}

/*!	\brief main function to call all the functionalities
 * 
 */

int main( int argc, char * argv[] )
{
	/*
	if( strcmp( argv[1] , "-Gauss_Curvature") == 0 && argc == 3 )
	{
		_Gauss_Curvature( argv[2] );
		return 0;
	}

	help( argv[0] );
	return 0;*/
	DelaunayTriangulation();
	Delaunay.write_m("DelaunayTriangulation");
	return 0;
}
 