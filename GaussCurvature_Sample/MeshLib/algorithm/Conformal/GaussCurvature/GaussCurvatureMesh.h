/*!
*      \file GaussCurvatureMesh.h
*      \brief Mesh for total Gaussian curvature
*	   \author David Gu
*      \date Documented 02/13/2014
*
*/

#ifndef  _GAUSS_CURVATURE_MESH_H_
#define  _GAUSS_CURVATURE_MESH_H_

#include <map>
#include <vector>

#include "Mesh/BaseMesh.h"
#include "Mesh/Vertex.h"
#include "Mesh/HalfEdge.h"
#include "Mesh/Edge.h"
#include "Mesh/Face.h"
#include "mesh/iterators.h"
#include "mesh/boundary.h"
#include "Parser/parser.h"

namespace MeshLib
{
/*!
*	\brief CGaussVertex class
*
*	Vertex class for Gauss Curvature
*/
class CGaussVertex : public CVertex
{

public:
	/*!
	 *	CHarmonicVertex constructor
	 */
	CGaussVertex() {m_k=0; };
	/*!
	 *	CHarmonicVertex destructor
	 */
	~CGaussVertex() {};
	/*
	 *	Gauss curvature
	 */
	double & k() { return m_k; };
	/*
	 *	Vertex normal
	 */
	CPoint & normal() { return m_normal; };

protected:
	/*
	 *	Gauss curvature
	 */
	double m_k;
	/*!
	 *	Vertex normal
	 */
	CPoint m_point;
};

/*!
*	\brief CGaussEdge class
*
*	Edge class for Gauss curvature
*/
class CGaussEdge : public  CEdge
  {
  public:
    /*!	CGaussEdge constructor
	 */
	 CGaussEdge() { m_length=0;  };
    /*!	CGaussEdge destructor
	 */
    ~CGaussEdge(){};
	/*! edge length trait
	 */
	double & length() { return m_length; };

  protected:
	/*! edge length trait */
	double   m_length;
};



/*!
*	\brief CGaussEdge class
*
*	Face class for Gauss curvature
*/
class CGaussFace : public  CFace
  {
  public:
    /*!	CGaussFace constructor
	 */
	 CGaussFace() {};
    /*!	CGaussFace destructor
	 */
    ~CGaussFace(){};
	/*!
	 *	face area
	 */
	double & area() { return m_area; };
	/*!
	 *	face normal
	 */
	CPoint & normal(){ return m_normal; };
 
  protected:
	double m_area;
	CPoint m_normal;
};


/*!
*	\brief CGaussHalfEdge class
*
*	HalfEdge class for Gauss curvature
*/
class CGaussHalfEdge : public  CHalfEdge
  {
  public:
    /*!	CHarmonicHalfEdge constructor
	 */
	 CGaussHalfEdge() {};
    /*!	CHarmonicHalfEdge destructor
	 */
    ~CGaussHalfEdge(){};
	/*!	Corner angle trait
	 */
	double & angle() { return m_angle; };

  protected:
	  /*! Corner angle trait */
	double m_angle;
};

/*-------------------------------------------------------------------------------------------------------------------------------------

	GaussCurvature Mesh Class

--------------------------------------------------------------------------------------------------------------------------------------*/
/*!
 *	\brief CGaussMesh class
 *
 *	Mesh class for Gauss curvature
 */
template<typename V, typename E, typename F, typename H>
class CGaussCurvatureMesh : public CBaseMesh<V,E,F,H>
{
public:
	
	typedef V CVertex;
	typedef E CEdge;
	typedef F CFace;
	typedef H CHalfEdge;

	typedef CGaussCurvatureMesh<V,E,F,H> M;

	typedef CBoundary<M> CBoundary;
	typedef CLoop<M> CLoop;
	
	typedef MeshVertexIterator<M> MeshVertexIterator;
	typedef MeshEdgeIterator<M> MeshEdgeIterator;
	typedef VertexVertexIterator<M> VertexVertexIterator;
	typedef VertexEdgeIterator<M> VertexEdgeIterator;
	typedef MeshFaceIterator<M> MeshFaceIterator;
	typedef FaceVertexIterator<M> FaceVertexIterator;
	typedef VertexFaceIterator<M> VertexFaceIterator;
	typedef FaceHalfedgeIterator<M> FaceHalfedgeIterator;
	typedef VertexInHalfedgeIterator<M> VertexInHalfedgeIterator;
};

/*! Mesh class for CHarmonicMapper class, Abbreviated as 'CGCMesh'
 */
typedef CGaussCurvatureMesh<CGaussVertex, CGaussEdge, CGaussFace, CGaussHalfEdge> CGCMesh;	

};
#endif  _GAUSS_CURVATURE_MESH_H_