/*!
*      \file Iterators.h
*      \brief Iterators for accessing geometric objects on a mesh
*	   \author David Gu
*      \date 10/07/2010
*
*/
 
#ifndef  _ITERATORS_H_
#define  _ITERATORS_H_

#include "BaseMesh.h"

namespace MeshLib{


//v->out halfedge
/*!
	\brief VertexOutHalfedgeIterator, transverse all the outgoing halfedges of a vertex ccwly.
*/
template<typename M>
class VertexOutHalfedgeIterator
{
public:
	/*!
	VertexOutHalfedgeIteartor constructor
	\param pMesh pointer to the current mesh
	\param v     pointer to the current vertex
	*/
	VertexOutHalfedgeIterator( M *  pMesh, typename M::CVertex *  v )
	{ m_pMesh = pMesh; m_vertex = v; m_halfedge = m_pMesh->vertexMostClwOutHalfEdge(v); };

	/*!
	VertexOutHalfedgeIterator destructor
	*/
	~VertexOutHalfedgeIterator(){};
	/*!
	prefix ++ operator, goes to the next ccw vertex out halfedge
	*/
	void operator++() //prefix
	{assert( m_halfedge != NULL ); 
	 if( m_halfedge == m_pMesh->vertexMostCcwOutHalfEdge(m_vertex) ) 
		 m_halfedge = NULL;
	 else
	 	 m_halfedge = m_pMesh->vertexNextCcwOutHalfEdge(m_halfedge); };

	/*!
		postfix ++ operator, goes to the next ccw vertex out halfedge
	*/
	void operator++(int) //postfix
	{assert( m_halfedge != NULL ); 
	 if( m_halfedge == m_pMesh->vertexMostCcwOutHalfEdge(m_vertex) ) 
		 m_halfedge = NULL;
	 else
	 	 m_halfedge = m_pMesh->vertexNextCcwOutHalfEdge(m_halfedge); };

	/*!
		The current halfedge the iterator pointing to.
	*/

	 typename M::CHalfEdge * value() { return m_halfedge; };
	 /*!
		whether all the out halfedges have been visited.
	 */
	 bool end(){ return m_halfedge == NULL; };
	/*!
		The current halfedge the iterator pointing to.
	*/
	 typename M::CHalfEdge * operator*() { return value(); };

private:
	/*!	
		Current mesh.
	*/
	M *        m_pMesh;
	/*!
		Current vertex.
	*/
	typename M::CVertex *      m_vertex;
	/*!
		Current halfedge.
	*/
	typename M::CHalfEdge * m_halfedge;
};

//v->in halfedge
/*!
	\brief VertexInHalfedgeIterator, transverse all the incoming halfedges of a vertex ccwly.
*/
template<typename M>
class VertexInHalfedgeIterator
{
public:
	/*!
	VertexInHalfedgeIteartor constructor
	\param pMesh pointer to the current mesh
	\param v     pointer to the current vertex
	*/
	VertexInHalfedgeIterator( M *  pMesh, typename M::CVertex * v )
	{ m_pMesh = pMesh; m_vertex = v; m_halfedge = m_pMesh->vertexMostClwInHalfEdge(v); };
	/*!
	VertexInHalfedgeIterator destructor
	*/
	~VertexInHalfedgeIterator(){};

	/*!
	prefix ++ operator, goes to the next ccw vertex in halfedge
	*/
		void operator++()	//prefix
	{
	 assert( m_halfedge != NULL ); 
	 if( m_halfedge == m_pMesh->vertexMostCcwInHalfEdge(m_vertex) ) 
		 m_halfedge = NULL; 
	 else
		m_halfedge = m_pMesh->vertexNextCcwInHalfEdge(m_halfedge); 
	};
	/*!
	postfix ++ operator, goes to the next ccw vertex in halfedge
	*/
	
	void operator++(int)	//postfix
	{
	 assert( m_halfedge != NULL ); 
	 if( m_halfedge == m_pMesh->vertexMostCcwInHalfEdge(m_vertex) ) 
		 m_halfedge = NULL; 
	 else
		m_halfedge = m_pMesh->vertexNextCcwInHalfEdge(m_halfedge); 
	};
	/*!
		The current halfedge the iterator pointing to.
	*/
	
	typename M::CHalfEdge * value() { return m_halfedge; };
	/*!	
		Indicate whether all the in halfedges of the vertex have been transversed.
	*/
	 bool end(){ return m_halfedge == NULL; };
	 /*!
		The current halfedge the iterator pointing to.
	 */
	 typename M::CHalfEdge * operator*() { return value(); };

private:
	/*!
		Current mesh.
	*/
	M *        m_pMesh;
	/*!
		Current vertex.
	*/
	typename M::CVertex *      m_vertex;
	/*!
		Current halfedge.
	*/
	typename M::CHalfEdge * m_halfedge;
};


//v -> v

/*!
	\brief VertexVertexIterator, transverse all the neighboring vertices of a vertex ccwly.
*/
template<typename M>
class VertexVertexIterator
{
public:
	/*!
		VertexVertexIterator constructor
		\param v the current vertex
	*/
	VertexVertexIterator( typename M::CVertex *  v )
	{ 
		m_vertex = v; 
		m_halfedge = (CHalfEdge*)m_vertex->most_clw_out_halfedge();
	};

	/*!
		VertexVertexIterator destructor
	*/
	~VertexVertexIterator(){};

	/*!
		VertexVertexIterator prefix operator ++, goes to the next neighboring vertex CCWly
	*/

	void operator++() //prefix
	{
		assert( m_halfedge != NULL ); 
		
		if( !m_vertex->boundary() )
		{
			if( m_halfedge != m_vertex->most_ccw_out_halfedge() )
			{
				m_halfedge = (CHalfEdge*)m_halfedge->ccw_rotate_about_source();
			}
			else
			{
				m_halfedge = NULL;
			}
			return;
		}

		if( m_vertex->boundary() )
		{
			if( m_halfedge == (CHalfEdge*)m_vertex->most_ccw_in_halfedge() )
			{
				m_halfedge = NULL;
				return;
			}

			CHalfEdge * he = (CHalfEdge*)m_halfedge->ccw_rotate_about_source();

			if( he == NULL )
			{
				m_halfedge = (CHalfEdge*)m_vertex->most_ccw_in_halfedge();
			}
			else
			{
				m_halfedge = he;
			}
		}

		return;
	};


	/*!
		VertexVertexIterator postfix operator ++, goes to the next neighboring vertex CCWly
	*/
	void operator++(int) //postfix
	{
		assert( m_halfedge != NULL ); 
		
		if( !m_vertex->boundary() )
		{
			if( m_halfedge != (CHalfEdge*)m_vertex->most_ccw_out_halfedge() )
			{
				m_halfedge = (CHalfEdge*)m_halfedge->ccw_rotate_about_source();
			}
			else
			{
				m_halfedge = NULL;
			}
			return;
		}

		if( m_vertex->boundary() )
		{
			if( m_halfedge == (CHalfEdge*)m_vertex->most_ccw_in_halfedge() )
			{
				m_halfedge = NULL;
				return;
			}

			CHalfEdge * he = (CHalfEdge*)m_halfedge->ccw_rotate_about_source();

			if( he == NULL )
			{
				m_halfedge = (CHalfEdge*)m_vertex->most_ccw_in_halfedge();
			}
			else
			{
				m_halfedge = he;
			}
		}

		return;
	};

	/*!
		The neighboring vertex, pointed by the current iterator
	*/

	typename M::CVertex * value() 
	 { 
		 if( !m_vertex->boundary() )
		 {
			 return (CVertex*)m_halfedge->target(); 
		 }

		 if( m_halfedge != (CHalfEdge*)m_vertex->most_ccw_in_halfedge() )
		 {
			 return (CVertex*)m_halfedge->target();
		 }

		 if( m_halfedge == (CHalfEdge*)m_vertex->most_ccw_in_halfedge() )
		 {
			 return (CVertex*)m_halfedge->source();
		 }
		 return NULL;
	 };

	/*!
		The neighboring vertex, pointed by the current iterator
	*/
	typename M::CVertex * operator*() { return value(); };

	/*!
		Indicate whether all the neighboring vertices have been accessed.
	*/
	 bool end(){ return m_halfedge == NULL; };

	/*!
		Reset the iterator.
	*/
	 void reset()	{ m_halfedge = (CHalfEdge*)m_vertex->most_clw_out_halfedge(); };

private:
	/*!
		Current vertex
	*/
	typename M::CVertex *   m_vertex;
	/*!	
		Current halfedge.
	*/
	typename M::CHalfEdge * m_halfedge;
};


//v -> edge
/*!
	\brief VertexEdgeIterator, transverse all the neighboring edges of a vertex ccwly.
*/

template<typename M>
class VertexEdgeIterator
{
public:
	/*!
		VertexEdgeIterator constructor
		\param v the current vertex
	*/
	VertexEdgeIterator( typename M::CVertex *  v )
	{ 
		m_vertex = v; 
		m_halfedge = (CHalfEdge*)m_vertex->most_clw_out_halfedge();
	};

	/*!
		VertexVertexIterator destructor
	*/
	~VertexEdgeIterator(){};
/*!
		VertexVertexIterator prefix operator ++, goes to the next neighboring vertex CCWly
	*/

	void operator++() //prefix
	{
		assert( m_halfedge != NULL ); 
		
		if( !m_vertex->boundary() )
		{
			if( m_halfedge != (CHalfEdge*)m_vertex->most_ccw_out_halfedge() )
			{
				m_halfedge = (CHalfEdge*)m_halfedge->ccw_rotate_about_source();
			}
			else
			{
				m_halfedge = NULL;
			}
			return;
		}

		if( m_vertex->boundary() )
		{
			if( m_halfedge == (CHalfEdge*)m_vertex->most_ccw_in_halfedge() )
			{
				m_halfedge = NULL;
				return;
			}

			M::CHalfEdge * he = (CHalfEdge*)m_halfedge->ccw_rotate_about_source();

			if( he == NULL )
			{
				m_halfedge = (CHalfEdge*)m_vertex->most_ccw_in_halfedge();
			}
			else
			{
				m_halfedge = he;
			}
		}

		return;
	};
	/*!
		VertexVertexIterator postfix operator ++, goes to the next neighboring vertex CCWly
	*/
	void operator++(int) //postfix
	{
		assert( m_halfedge != NULL ); 
		
		if( !m_vertex->boundary() )
		{
			if( m_halfedge != (CHalfEdge*)m_vertex->most_ccw_out_halfedge() )
			{
				m_halfedge = (CHalfEdge*)m_halfedge->ccw_rotate_about_source();
			}
			else
			{
				m_halfedge = NULL;
			}
			return;
		}

		if( m_vertex->boundary() )
		{
			if( m_halfedge == (CHalfEdge*)m_vertex->most_ccw_in_halfedge() )
			{
				m_halfedge = NULL;
				return;
			}

			M::CHalfEdge * he = (CHalfEdge*)m_halfedge->ccw_rotate_about_source();

			if( he == NULL )
			{
				m_halfedge = (CHalfEdge*)m_vertex->most_ccw_in_halfedge();
			}
			else
			{
				m_halfedge = he;
			}
		}

		return;
	};
	/*!
		The neighboring edge, pointed by the current iterator
	*/

	typename M::CEdge * value() 
	 { 
		 assert( m_halfedge != NULL );
		 return (CEdge*)m_halfedge->edge();
	 };

	/*!
		The neighboring edge, pointed by the current iterator
	*/
	typename M::CEdge * operator*() { return value(); };
	/*!
		Indicate whether all the neighboring edges have been accessed.
	*/
	 bool end(){ return m_halfedge == NULL; };
	/*!
		Reset the VerexEdgeIterator.
	*/
	 void reset()	{ m_halfedge = (CHalfEdge*)m_vertex->most_clw_out_halfedge(); };

private:
	/*! current vertex 
	*/
	typename M::CVertex *   m_vertex;
	/*! current halfedge
	*/
	typename M::CHalfEdge * m_halfedge;
};



// v->face
/*!
	\brief VertexFaceIterator, transverse all the neighboring faces of a vertex ccwly.
*/
template<typename M>
class VertexFaceIterator
{
public:
	/*!
		VertexFaceIterator constructor
		\param v the current vertex
	*/
	VertexFaceIterator( typename M::CVertex * & v )
	{ 
		m_vertex = v; 
		m_halfedge = (M::CHalfEdge*)m_vertex->most_clw_out_halfedge(); 
	};
	/*!
		VertexFaceIterator destructor
	*/
	~VertexFaceIterator(){};
	/*!
		VertexFaceIterator prefix operator ++, goes to the next neighboring face CCWly
	*/
	void operator++() //prefix
	{
		assert( m_halfedge != NULL );  

		if( m_halfedge == (M::CHalfEdge*)m_vertex->most_ccw_out_halfedge() ) 
		{
			m_halfedge = NULL;
			return;
		}
		m_halfedge = (M::CHalfEdge*)m_halfedge->ccw_rotate_about_source();
	};
	/*!
		VertexFaceIterator prefix operator ++, goes to the next neighboring face CCWly
	*/

	void operator++(int) //postfix
	{
		assert( m_halfedge != NULL );  

		if( m_halfedge == (M::CHalfEdge*)m_vertex->most_ccw_out_halfedge() ) 
		{
			m_halfedge = NULL;
			return;
		}
		m_halfedge = (M::CHalfEdge*)m_halfedge->ccw_rotate_about_source();
	};
	/*!
		The neighboring face, pointed by the current iterator
	*/
	typename M::CFace * value() { return (M::CFace*) m_halfedge->face(); };
	/*!
		The neighboring face, pointed by the current iterator
	*/
	typename M::CFace * operator*() { return value(); };
	/*!
		Indicate whether all the neighboring faces have been accessed.
	*/
	 bool end(){ return m_halfedge == NULL; };
	 /*!
	 Reset the VertexFaceIterator
	 */
	 void reset()	{ m_halfedge = (M::CHalfEdge*)m_vertex->most_clw_out_halfedge(); };

private:
	/*!
	current vertex
	*/
	typename M::CVertex *   m_vertex;
	/*!
	current halfedge
	*/
	typename M::CHalfEdge * m_halfedge;
};

// f -> halfedge
/*!
	\brief FaceHalfedgeIterator, transverse all the halfedges of a face CCWly.
*/

template<typename M>
class FaceHalfedgeIterator
{
public:
	/*!
		FaceHalfedgeIterator constructor
		\param f the current face
	*/
	FaceHalfedgeIterator( typename M::CFace * f )
	{ 
		m_face = f; 
		m_halfedge = (M::CHalfEdge*)f->halfedge(); 
	};
	/*!
		FaceHalfedgeIterator destructor
	*/
	~FaceHalfedgeIterator(){};
	/*!
		VertexVertexIterator prefix operator ++, goes to the next halfedge CCWly
	*/
	void operator++() //prefix
	{
		assert( m_halfedge != NULL );
		m_halfedge = (M::CHalfEdge*)m_halfedge->he_next();

		if( m_halfedge == m_face->halfedge() )
		{
			 m_halfedge = NULL;
			return;
		};
	}

	/*!
		VertexVertexIterator prefix operator ++, goes to the next halfedge CCWly
	*/
	void operator++(int) //postfix
	{
		assert( m_halfedge != NULL );
		m_halfedge = (M::CHalfEdge*)m_halfedge->he_next();

		if( m_halfedge == m_face->halfedge() )
		{
			 m_halfedge = NULL;
			return;
		};
	}

	/*!
		The halfedge, pointed by the current iterator
	*/
	typename M::CHalfEdge * value() { return m_halfedge; };
	/*!
		The halfedge, pointed by the current iterator
	*/
	typename M::CHalfEdge * operator*() { return value(); };

	/*!
		Indicate whether all the halfedges have been accessed.
	*/
	bool end(){ return m_halfedge == NULL; };

private:
	/*!
		current face
	*/
	typename M::CFace *        m_face;
	/*!
		current halfedge
	*/
	typename M::CHalfEdge * m_halfedge;
};


// f -> edge
/*!
	\brief FaceEdgeIterator, transverse all the edges of a face CCWly.
*/
template<typename M>
class FaceEdgeIterator
{
public:
	/*!
		FaceEdgeIterator constructor
		\param f the current face
	*/	
	FaceEdgeIterator( typename M::CFace * f )
	{ 
		m_face = f; 
		m_halfedge = (M::CHalfEdge*)f->halfedge(); 
	};

	/*!
		FaceEdgeIterator destructor
	*/	
	~FaceEdgeIterator(){};
	/*!
		FaceEdgeIterator prefix operator ++, goes to the next edge CCWly
	*/
	void operator++()	//prefix
	{
		assert( m_halfedge != NULL );
		m_halfedge = (M::CHalfEdge*)m_halfedge->he_next();

		if( m_halfedge == (M::CHalfEdge*)m_face->halfedge() )
		{
			 m_halfedge = NULL;
			return;
		};
	}

	/*!
		FaceEdgeIterator prefix operator ++, goes to the next edge CCWly
	*/
	void operator++(int)	//postfix
	{
		assert( m_halfedge != NULL );
		m_halfedge = (M::CHalfEdge*)m_halfedge->he_next();

		if( m_halfedge == m_face->halfedge() )
		{
			 m_halfedge = NULL;
			return;
		};
	}
	/*!
		The edge, pointed by the current iterator
	*/
	typename M::CEdge * value() { return (M::CEdge*) m_halfedge->edge(); };
	/*!
		The edge, pointed by the current iterator
	*/
	typename M::CEdge * operator*() { return value(); };
	/*!
		Indicate whether all the edges have been transversed.
	*/
	bool end(){ return m_halfedge == NULL; };

private:
	/*! Current face. */
	typename M::CFace  *       m_face;
	/*! Current halfedge. */
	typename M::CHalfEdge * m_halfedge;
};


// f -> vertex
/*!
	\brief FaceVertexIterator, transverse all the vertices of a face CCWly.
*/
template<typename M>
class FaceVertexIterator
{
public:
	/*!
		FaceVertexIterator constructor
		\param f the current face
	*/
	FaceVertexIterator( typename M::CFace * f )
	{ 
		m_face = f; 
		m_halfedge = (M::CHalfEdge*)f->halfedge(); 
	};
	/*!
		FaceVertexIterator destructor
	*/
	
	~FaceVertexIterator(){};
	/*!
		FaceVertexIterator prefix operator ++, goes to the next vertex CCWly
	*/
	void operator++() //prefix
	{
		assert( m_halfedge != NULL );
		m_halfedge = (M::CHalfEdge*)m_halfedge->he_next();

		if( m_halfedge == (M::CHalfEdge*)m_face->halfedge() )
		{
			 m_halfedge = NULL;
			return;
		};
	}

	/*!
		FaceVertexIterator prefix operator ++, goes to the next vertex CCWly
	*/
	void operator++(int) //postfix
	{
		assert( m_halfedge != NULL );
		m_halfedge = (M::CHalfEdge*)m_halfedge->he_next();

		if( m_halfedge == (M::CHalfEdge*)m_face->halfedge() )
		{
			 m_halfedge = NULL;
			return;
		};
	}
	/*!
		The vertex, pointed by the current iterator
	*/
	typename M::CVertex * value() { return (M::CVertex*) m_halfedge->target(); };
	/*!
		The vertex, pointed by the current iterator
	*/
	typename M::CVertex * operator*() { return value(); };
	/*!
		Indicate whether all the vertices have been accessed.
	*/
	bool end(){ return m_halfedge == NULL; };

private:
	/*!	Current face.
	*/
	typename M::CFace         * m_face;
	/*!	Current halfedge.
	*/
	typename M::CHalfEdge * m_halfedge;
};


// mesh->v
/*!
	\brief MeshVertexIterator, transverse all the vertices in the mesh.
*/

template<typename M>
class MeshVertexIterator
{
public:
	/*!
	MeshVertexIterator constructor, 
	\param pMesh the current mesh
	*/
	MeshVertexIterator( M * pMesh )
	{
		m_pMesh = pMesh;
		m_iter = m_pMesh->vertices().begin();
	}
	/*!
	The vertex, pointed by the current iterator
	*/
	typename M::CVertex * value() { return *m_iter; };
	/*!
	The vertex, pointed by the current iterator
	*/
		
	typename M::CVertex * operator*(){ return value(); };
	/*!
		MeshVertexIterator prefix operator ++, goes to the next vertex 
	*/	
	void operator++()	 { ++ m_iter; }; //prefix
	/*!
		MeshVertexIterator prefix operator ++, goes to the next vertex 
	*/	
	void operator++(int) { ++ m_iter; }; //postfix
	/*!
		Indicate whether all the vertices have been accessed.
	*/
	bool end() { return m_iter == m_pMesh->vertices().end(); }
	
private:
	/*!
		Current mesh.
	*/
	M * m_pMesh;
	/*! 
	Current vertex list iterator.
	*/
	typename std::list<typename M::CVertex*>::iterator m_iter;
};

// mesh->f
/*!
	\brief MeshFaceIterator, transverse all the faces in the mesh.
*/
template<typename M>
class MeshFaceIterator
{
public:
	/*!
	MeshFaceIterator constructor, 
	\param pMesh the current mesh
	*/
	MeshFaceIterator( M * pMesh )
	{
      m_pMesh = pMesh;
      m_iter = pMesh->faces().begin();
	}
	/*!
	The face, pointed by the current iterator
	*/
	typename M::CFace * value() { return *m_iter; };
	/*!
	The face, pointed by the current iterator
	*/
	typename M::CFace * operator*(){ return value(); };

	/*!
		MeshFaceIterator prefix operator ++, goes to the next vertex 
	*/
	void operator++() { ++ m_iter; }; //prefix
	/*!
		MeshFaceIterator postfix operator ++, goes to the next vertex 
	*/
	void operator++(int) { ++ m_iter; }; //postfix
	/*!
		Indicate whether all the faces have been accessed.
	*/
	bool end() { return m_iter == m_pMesh->faces().end(); }

private:
	/*! Current mesh.
	*/
	M * m_pMesh;
	/*! Current face list iterator.
	*/
	typename std::list<typename M::CFace*>::iterator  m_iter;
};

//Mesh->e
/*!
	\brief MeshEdgeIterator, transverse all the edges in the mesh.
*/
template<typename M>
class MeshEdgeIterator
{
public:
	/*!
	MeshEdgeIterator constructor, 
	\param pMesh the current mesh
	*/	
	MeshEdgeIterator( M * pMesh )
	{
		m_pMesh = pMesh;
		m_iter = m_pMesh->edges().begin();
	}
	/*!
	The edge, pointed by the current iterator
	*/	
	typename M::CEdge * value() { return *m_iter; };
	/*!
	The edge, pointed by the current iterator
	*/	
	typename M::CEdge * operator*(){ return value(); };
	/*!
		MeshEdgeIterator prefix operator ++, goes to the next edge
	*/	
	void operator++() { ++ m_iter; }; //prefix
	/*!
		MeshEdgeIterator postfix operator ++, goes to the next edge
	*/	
	void operator++(int) {m_iter++; }; //postfix
	/*!
		Indicate whether all the edges have been accessed.
	*/	
	bool end() { return m_iter == m_pMesh->edges().end(); }


private:
	/*!
	current mesh
	*/
	M * m_pMesh;
	/*!
	current edge list iterator
	*/
	typename std::list<typename M::CEdge*>::iterator m_iter;
};

// Mesh->he
/*!
	\brief MeshHalfEdgeIterator, transverse all the halfedges in the mesh.
*/
template<typename M>
class MeshHalfEdgeIterator
{
public:
	/*!
	MeshHalfedgeIterator constructor, 
	\param pMesh the current mesh
	*/
	MeshHalfEdgeIterator( M * pMesh )
	{
     m_pMesh = pMesh;
     m_iter = m_pMesh->edges().begin();
     m_id = 0;
	}
	/*!
	The halfedge, pointed by the current iterator
	*/	
	typename M::CHalfEdge * value() { typename M::CEdge * e = *m_iter; return (typename M::CHalfEdge*)e->halfedge(m_id); };
	/*!
	The halfedge, pointed by the current iterator
	*/	
	typename M::CHalfEdge * operator*(){ return value(); };
	/*!
		MeshVertexIterator prefix operator ++, goes to the next halfedge 
	*/
	void operator++() //prefix
	{ 
		++m_id;

		switch( m_id )
		{
		case 1:
			{
				M::CEdge * e = *m_iter;
				if( e->halfedge(m_id) == NULL )
				{
					m_id = 0;
					++ m_iter;
				}
			}
			break;
		case 2:
			m_id = 0;
			++m_iter;
			break;
		}
	};
	/*!
		MeshVertexIterator postfix operator ++, goes to the next vertex 
	*/
	void operator++(int) //postfix
	{ 
		++m_id;

		switch( m_id )
		{
		case 1:
			{
				M::CEdge * e = *m_iter;
				if( e->halfedge(m_id) == NULL )
				{
					m_id = 0;
					++ m_iter;
				}
			}
			break;
		case 2:
			m_id = 0;
			++m_iter;
			break;
		}
	};
	/*!
	Indicate whether all the halfedges have been accessed
	*/
	bool end() { return m_iter == m_pMesh->edges().end(); }
	

private:
	/*!
		Current halfedge
	*/
	typename M::CHalfEdge * m_he;
	/*!
		Current mesh
	*/
	M *	 m_pMesh;
	/*!
		Current edge list iterator
	*/
	typename std::list<typename M::CEdge*>::iterator m_iter;
	int  m_id;
};


} //Iterators

#endif
