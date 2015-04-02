/*!
*      \file Boundary.h
*      \brief Trace boundary loops
*	   \author David Gu
*      \date 10/07/2010
*
*/


#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <list>
#include <set>

#include "../Mesh/BaseMesh.h"
#include "../Mesh/iterators.h"

namespace MeshLib
{
/*!
	\brief CLoop Boundary loop  class.	
*/

template<typename M>
class CLoop
{
public:
	/*!
		Constructor of the CLoop
		\param pMesh  pointer to the current mesh
		\param pH halfedge on the boundary loop
	*/
	CLoop( M * pMesh, typename M::CHalfEdge * pH );
	/*!
		Constructor of the CLoop
		\param pMesh  pointer to the current mesh
	*/
	 CLoop( M * pMesh ) { m_pMesh = pMesh; m_length = 0; m_pHalfedge = NULL; };
	 /*!
		Destructor of CLoop.
	 */
	~CLoop();
	
	/*!
	The list of haledges on the current boundary loop.
	*/
	std::list<typename M::CHalfEdge*>  & halfedges() 
	{
		return m_halfedges;
	}
	
	/*!
		The length of the current boundary loop.
	*/
	double length() 
	{ 
		return m_length; 
	}
	/*!
	 *  Output to a file
	 */
	void write( const char * file );
	/*!
	 *  Input from a file
	 */
	void read( const char * file );

protected:
	/*!
		Pointer to the current mesh.
	*/
	M		* m_pMesh;
	/*! The length of the current boundary loop.
	*/
	double										    m_length;
	/*!
		Starting halfedge of the current boundary loop.
	*/
	typename M::CHalfEdge						  * m_pHalfedge;
	/*!
		List of consecutive halfedges along the boundary loop.
	*/
	std::list<typename M::CHalfEdge*>				m_halfedges;

};

/*!
	\brief CBoundary Boundary  class.	
	\tparam CVertex Vertex type
	\tparam CEdge   Edge   type
	\tparam CFace   Face   type
	\tparam CHalfEdge HalfEdge type
*/

template<typename M>
class CBoundary
{

public:
	/*!
	CBoundary constructor
	\param pMesh pointer to the current mesh
	*/
	CBoundary( M * pMesh );
	/*!
	CBoundary destructor
	*/
	~CBoundary();
	/*!
	The list of boundary loops.
	*/
	std::vector<CLoop<M>*> & loops() 
	{ 
		return m_loops; 
	} 

protected:
	/*!
		Pointer to the current mesh.
	*/
	M * m_pMesh;
	/*!
		List of boundary loops.
	*/
	typename std::vector<CLoop<M>*> m_loops;
	/*!
		Bubble sort the loops
		\param loops the vector of loops
	*/
	void _bubble_sort( std::vector<CLoop<M>*> & loops);
};

/*!
	CLoop constructure, trace the boundary loop, starting from the halfedge pH.
	\param pMesh pointer to the current mesh
	\param pH  halfedge on the current boundary loop
*/
template<typename M>
CLoop<M>::CLoop( M *pMesh, typename M::CHalfEdge * pH )
{
	m_pMesh     = pMesh;
	m_pHalfedge = pH;

	m_length = 0;
	M::CHalfEdge * he = pH;

	//trace the boundary loop
	do{
		M::CVertex * v = m_pMesh->halfedgeTarget( he ); //(M::CVertex*)he->target();
		he = m_pMesh->vertexMostClwOutHalfEdge( v );
		m_halfedges.push_back( he );
		m_length += m_pMesh->edgeLength( (M::CEdge*)he->edge() );
	}while( he != m_pHalfedge );
}

/*!
CLoop destructor, clean up the list of halfedges in the loop
*/
template<typename M>
CLoop<M>::~CLoop()
{
	m_halfedges.clear();
}


/*!
	_bubble_sort
	bubble sort  a vector of boundary loop objects, according to their lengths
	\param loops vector of loops
*/
template<typename M>
void CBoundary<M>::_bubble_sort( std::vector<CLoop<M>*> & loops)
{
      int i, j, flag = 1;    // set flag to 1 to start first pass
      CLoop<M> * temp;             // holding variable
      int numLength = (int)loops.size( ); 
      for(i = 1; (i <= numLength) && flag; i++)
     {
          flag = 0;
          for (j=0; j < (numLength -1); j++)
         {
               if (loops[j+1]->length() > loops[j]->length() )      // ascending order simply changes to <
              { 
                    temp = loops[j];								// swap elements
                    loops[j] = loops[j+1];
                    loops[j+1] = temp;
                    flag = 1;										// indicates that a swap occurred.
               }
          }
     }
}

/*!
	CBoundary constructor
	\param pMesh the current mesh
*/
template<typename M>
CBoundary<M>::CBoundary( M * pMesh )
{
	m_pMesh = pMesh;
	//collect all boundary halfedges
	std::set<M::CHalfEdge*> boundary_hes;
	for( M::MeshEdgeIterator eiter( m_pMesh); !eiter.end(); eiter ++ )
	{
		M::CEdge * e = *eiter;
		if( !m_pMesh->isBoundary(e) ) continue;

		M::CHalfEdge * he = m_pMesh->edgeHalfedge( e, 0);
		boundary_hes.insert( he );
	}
	//trace all the boundary loops
	while( !boundary_hes.empty() )
	{
		//get the first boundary halfedge
		std::set<M::CHalfEdge*>::iterator siter = boundary_hes.begin();
		M::CHalfEdge * he = *siter;
		//trace along this boundary halfedge
		CLoop<M> * pL = new CLoop<M>( m_pMesh, he );
		assert(pL);
		m_loops.push_back( pL );
		//remove all the boundary halfedges, which are in the same boundary loop as the head, from the halfedge list
		for( std::list<M::CHalfEdge*>::iterator hiter = pL->halfedges().begin(); hiter != pL->halfedges().end(); hiter ++ )
		{
			M::CHalfEdge * he = *hiter;
			siter = boundary_hes.find( he );
			if( siter == boundary_hes.end() ) continue;
			boundary_hes.erase( siter );
		}
	}
	
	//std::sort( vlps.begin(), vlps.end(), loop_compare<CVertex,CFace,CEdge,CHalfEdge> );
	_bubble_sort( m_loops );
}

/*!	CBoundary destructor, delete all boundary loop objects.
*/
template<typename M>
CBoundary<M>::~CBoundary()
{
	for( int i = 0; i < (int)m_loops.size(); i ++ )
	{
		CLoop<M> * pL = m_loops[i];
		delete pL;
	}
};

/*!
	Output the loop to a file
	\param file_name the name of the file
*/
template<typename M>
void CLoop<M>::write( const char * file_name )
{
	std::ofstream myfile;
	myfile.open (file_name);
	for( std::list<M::CHalfEdge*>::iterator hiter = m_halfedges.begin(); hiter != m_halfedges.end(); hiter ++ )
	{
		M::CHalfEdge * pH = *hiter;
		M::CVertex * pV = m_pMesh->halfedgeSource(pH);
		M::CVertex * pW = m_pMesh->halfedgeTarget(pH);
		myfile << pV->id() << " " << pW->id() << std::endl;
	}
	myfile.close();
};

/*!
	Output the loop to a file
	\param file_name the name of the file
*/
template<typename M>
void CLoop<M>::read( const char * file_name )
{
	std::ifstream myfile;
	myfile.open (file_name);

	if (myfile.is_open())
	{
		while ( myfile.good() )
		{
			int source, target;
			myfile >> source >> target;
			M::CVertex * pS = m_pMesh->idVertex( source );
			M::CVertex * pT = m_pMesh->idVertex( target );
			M::CEdge   * pE = m_pMesh->vertexEdge( pS, pT );
			M::CHalfEdge * pH = m_pMesh->edgeHalfedge(pE,0);
			m_halfedges.push_back( pH );
		}
		myfile.close();
	}

};


}
#endif

