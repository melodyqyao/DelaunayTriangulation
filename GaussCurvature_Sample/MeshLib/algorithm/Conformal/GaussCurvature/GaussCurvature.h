/*!
*      \file GaussCurvature.h
*      \brief Algorithm for Gauss Curvature
*	   \author David Gu
*      \date Document 02/13/2014
*
*/


#ifndef _GAUSS_CURVATURE_H_
#define _GAUSS_CURVATURE_H_

#include <vector>
#include "Mesh/iterators.h"
#include "GaussCurvatureMesh.h"

#ifndef PI
#define PI 3.14159265358979323846
#endif
namespace MeshLib
{

	template<typename M>
	class CGaussCurvature
	{
	public:
		/*!	CGaussCurvature constructor
		 *	\param pMesh the input mesh
		 */
		CGaussCurvature( M* pMesh);
		/*!	CHarmonicMapper destructor
		 */
		~CGaussCurvature();
		/*!  Compute the harmonic map using direct method
		 */
		void _calculate_curvature();
		/*!	Compute face normal
		 */
		void _calculate_face_normal();
		/*!	Compute vertex normal
		 */
		void _calculate_vertex_normal();
		/*! Compute Euler number
		 */
		void _calculate_Euler_characteristics();

	protected:
		/*!	The input surface mesh
		 */
		M* m_pMesh;
		/*	boundary
		 */
		typename M::CBoundary m_boundary;
	};

double _cosine_law( double a, double b, double c )
{
          double cs =  ( a*a + b * b  - c * c )/( 2.0 * a * b );
          assert( cs <= 1.0 && cs >= -1.0 );
          return acos( cs );    
};


/*!	CHarmonicMapper constructor 
*	Count the number of interior vertices, boundary vertices and the edge weight
*
*/
template<typename M>
CGaussCurvature<M>::CGaussCurvature( M* pMesh ): m_pMesh( pMesh ), m_boundary( pMesh )
{
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
	{
		M::CVertex * pV = *viter;
		pV->k() = 0;
	}
};

/*!
 *	CGaussCurvature destructor
 */
template<typename M>
CGaussCurvature<M>::~CGaussCurvature()
{
};


//Set the boundary vertices to the unit circle
/*!
 *	Fix the boundary using arc length parameter
 */
template<typename M>
void CGaussCurvature<M>::_calculate_curvature()
{
	//insert your code here	
	double sum=0.0;
	double k_i;
	//calculate edge length
	for(M::MeshEdgeIterator eiter(m_pMesh);!eiter.end();++eiter)
	{
		M::CEdge *pE = *eiter;
		M::CVertex *pV1 = m_pMesh->edgeVertex1(pE);
		M::CVertex *pV2 = m_pMesh->edgeVertex2(pE);
		pE->length() = (pV1->point() - pV2->point()).norm();
	}
	//calculate corner angle
    //std::cout << "I am here!"<< std::endl;;
	for(M::MeshFaceIterator fiter(m_pMesh);!fiter.end();++fiter)
	{
			M::CFace *pF = *fiter;
			std::vector<M::CHalfEdge*> v;
			std::vector<double> len;
		for(M::FaceHalfedgeIterator fhiter(pF);!fhiter.end();++fhiter)
		{
			M::CHalfEdge *he = *fhiter;
			v.push_back(he); 
			len.push_back(((CGaussEdge*)he->edge())->length());
		}	
			 v[0]->angle() = _cosine_law(len[0],len[1],len[2]);
			 v[1]->angle() = _cosine_law(len[1],len[2],len[0]);
			 v[2]->angle() = _cosine_law(len[0],len[2],len[1]);	
			 v.clear();
			 len.clear();
			// std::cout<<"I am here!"<<std::endl;
	}

	//caculate curvature
	for(M::MeshVertexIterator viter(m_pMesh); !viter.end();++viter)
	{
		M::CVertex *pV = *viter;
		//if vertex is on the boundary
		if(pV->boundary()){
			 k_i = PI;
			for(M::VertexInHalfedgeIterator vihiter(m_pMesh,pV);!vihiter.end();++vihiter)
			{
				M::CHalfEdge *ihe = *vihiter;
				k_i -= ihe->angle();
			}
		}
		//if vertex is not on the boundary
		else{
			 k_i = 2*PI;
			for(M::VertexInHalfedgeIterator vihiter(m_pMesh,pV);!vihiter.end();++vihiter)
			{
				M::CHalfEdge *ihe = *vihiter;
				k_i -= ihe->angle();
			}
		}
		sum += k_i; 
	}
	std::cout << "Total Curvature is " <<  sum/PI << " PI" << std::endl;;
}

/*!
 *	Compute face normal
 */
template<typename M>
void CGaussCurvature<M>::_calculate_face_normal()
{
	//insert your code here
	for(M::MeshFaceIterator fiter(m_pMesh);!fiter.end();++fiter)
	{

		M::CFace *pF=*fiter;	
		CPoint d;
		std::vector<CPoint> v;
		for(M::FaceVertexIterator fviter(pF);!fviter.end();++fviter)
		{
			M::CVertex *pV=*fviter;		
			CPoint p= pV->point();
			v.push_back(p);
		}
			d=(v[1]-v[0])^(v[2]-v[0]);
			pF->area() = 1/2*d.norm();
			pF->normal() = d/d.norm();
			v.clear();
	}
}

/*!
 *	Compute vertex normal
 */
template<typename M>
void CGaussCurvature<M>::_calculate_vertex_normal()
{
	//insert your code here 
	for(M::MeshVertexIterator viter(m_pMesh);!viter.end();++viter)
	{
		M::CVertex *pV=*viter;
		CPoint d,d_i; //d is the sum of the normal, d_i is the average normal
		double s=0.0; //the area for each face
		for(M::VertexFaceIterator vfiter(pV);!vfiter.end();++vfiter)
		{
			M::CFace *pF= *vfiter;
			d += (pF->normal())*(pF->area());
			s += pF->area();
		}
		d_i = d/s;
		pV->normal()=d_i/d_i.norm();
	}
}

/*!
 *  Calculate Euler characteristics number
 */
template<typename M>
void CGaussCurvature<M>::_calculate_Euler_characteristics()
{
	//insert your code in this function
	int V, E, F;
	V = m_pMesh->numVertices();
	E = m_pMesh->numEdges();
	F = m_pMesh->numFaces();
	std::cout << "Vertices: " << V << " Faces: " << F << " Edges: " << E << std::endl;
	int euler = V + F - E;
	std::cout << "Euler Characteristic Number " << euler << std::endl;	
	int b;
	int g;
	b = M::CBoundary(m_boundary).loops().size(); 
	std::cout << "Number of boundaries " << b << std::endl;
	g = (2-b-euler)/2;
	std::cout << "Genus " << g << std::endl;
}
};
#endif

