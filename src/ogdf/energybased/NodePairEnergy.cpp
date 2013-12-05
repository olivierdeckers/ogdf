/*
 * $Revision: 2813 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2012-10-13 14:05:35 +0200 (Sa, 13. Okt 2012) $
 ***************************************************************/

/** \file
 * \brief Implementation of class NodePairEnergy
 *
 * \author Rene Weiskircher
 *
 * \par License:
 * This file is part of the Open Graph Drawing Framework (OGDF).
 *
 * \par
 * Copyright (C)<br>
 * See README.txt in the root directory of the OGDF installation for details.
 *
 * \par
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * Version 2 or 3 as published by the Free Software Foundation;
 * see the file LICENSE.txt included in the packaging of this file
 * for details.
 *
 * \par
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * \par
 * You should have received a copy of the GNU General Public
 * License along with this program; if not, write to the Free
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 *
 * \see  http://www.gnu.org/copyleft/gpl.html
 ***************************************************************/


#include <ogdf/internal/energybased/NodePairEnergy.h>


namespace ogdf {


NodePairEnergy::NodePairEnergy(const string energyname, GraphAttributes &AG) :
	EnergyFunction(energyname,AG),
	m_shape(m_G),
	m_adjacentOracle(m_G)
{
	node v;
	double lengthSum = 0;
	forall_nodes(v,m_G) { //saving the shapes of the nodes in m_shape
		DPoint center(AG.x(v),AG.y(v));
		lengthSum += AG.width(v);
		lengthSum += AG.height(v);
		m_shape[v] = IntersectionRectangle(center,AG.width(v),AG.height(v));
	}
	m_G.allNodes(m_nonIsolated);
	ListIterator<node> it, itSucc;
	for(it = m_nonIsolated.begin(); it.valid(); it = itSucc) {
		itSucc = it.succ();
		if((*it)->degree() == 0) m_nonIsolated.del(it);
	}
	m_nodeNums = OGDF_NEW NodeArray<int>(m_G,0);
	int n_num = 1;
	for(it = m_nonIsolated.begin(); it.valid(); ++it) {
		(*m_nodeNums)[*it] = n_num;
		n_num++;
	}
	n_num--;
	m_pairEnergy = new Array2D<double> (1,n_num,1,n_num);
	m_candPairEnergy = new Array2D<double> (1,n_num,1,n_num);
}


void NodePairEnergy::computeEnergy()
{
	int n_num = m_nonIsolated.size();
	double energySum = 0.0;
	Array<node> numNodes(1,n_num);

	ListIterator<node> it;
	for(it = m_nonIsolated.begin(); it.valid(); ++it) {
		numNodes[(*m_nodeNums)[*it]] = *it;
	}
	for(int i = 1; i <= n_num-1 ; i++) {
		for(int j = i+1; j <= n_num; j++) {
			double E = computePairEnergy(numNodes[i],numNodes[j]);
			(*m_pairEnergy)(i,j) = E;
			energySum += E;
		}
	}
	m_energy = energySum;
}


double NodePairEnergy::computePairEnergy(const node v, const node w) const {
	return computeCoordEnergy(v,w,currentPos(v),currentPos(w));
}


void NodePairEnergy::internalCandidateTaken() {
	node s = sourceTestNode();
	internalCandidateTaken(s);
	node t = targetTestNode();
	if(t != NULL) 
		internalCandidateTaken(t);
}

void NodePairEnergy::internalCandidateTaken(const node n) {
	int candNum = (*m_nodeNums)[n];
	ListIterator<node> it;
	for(it = m_nonIsolated.begin(); it.valid(); ++ it) {
		if((*it) != n) {
			int numit = (*m_nodeNums)[*it];
			(*m_pairEnergy)(min(numit,candNum),max(numit,candNum)) = (*m_candPairEnergy)(min(numit,candNum),max(numit,candNum));
			(*m_candPairEnergy)(min(numit,candNum),max(numit,candNum)) = 0.0;
		}
	}
}


void NodePairEnergy::compCandEnergy()
{
	node s = sourceTestNode();
	node t = targetTestNode();
	int nums = (*m_nodeNums)[s];
	int numt = (*m_nodeNums)[t];

	compCandEnergy(s, t, sourceTestPos());

	if(t != NULL) {
		compCandEnergy(t, s, targetTestPos());

		m_candidateEnergy -= (*m_pairEnergy)(min(nums,numt),max(nums,numt));
		double coordEnergy = computeCoordEnergy(s,t,sourceTestPos(),targetTestPos());
		(*m_candPairEnergy)(min(nums,numt),max(nums,numt)) = coordEnergy;
		m_candidateEnergy += coordEnergy;

		if(m_candidateEnergy < 0.0) {
			OGDF_ASSERT(m_candidateEnergy > -0.00001);
			m_candidateEnergy = 0.0;
		}
	}
	OGDF_ASSERT(m_candidateEnergy >= -0.0001);
}

void NodePairEnergy::compCandEnergy(const node v, const node ignore, const DPoint newPos)
{
	int numv = (*m_nodeNums)[v];
	m_candidateEnergy = energy();
	ListIterator<node> it;
	for(it = m_nonIsolated.begin(); it.valid(); ++ it) {
		int j = (*m_nodeNums)[*it];
		if(*it != v && *it != ignore) {
			m_candidateEnergy -= (*m_pairEnergy)(min(j,numv),max(j,numv));
			double coordEnergy = computeCoordEnergy(v,*it,newPos,currentPos(*it));
			(*m_candPairEnergy)(min(j,numv),max(j,numv)) = coordEnergy;
			m_candidateEnergy += coordEnergy;
			if(m_candidateEnergy < 0.0) {
				OGDF_ASSERT(m_candidateEnergy > -0.00001);
				m_candidateEnergy = 0.0;
			}
		}
		else (*m_candPairEnergy)(min(j,numv),max(j,numv)) = 0.0;
	}
}


#ifdef OGDF_DEBUG
void NodePairEnergy::printInternalData() const {
	cout << "\nCandidate energies:";
	for(int i=1; i< m_nonIsolated.size(); i++)
		for(int j=i+1; j <= m_nonIsolated.size(); j++)
			if((*m_candPairEnergy)(i,j) != 0.0)
				cout << "\nCandidateEnergy(" << i << ',' << j << ") = " << (*m_candPairEnergy)(i,j);
	cout << "\nPair energies:";
	for(int i=1; i< m_nonIsolated.size(); i++)
		for(int j=i+1; j <= m_nonIsolated.size(); j++)
			if((*m_pairEnergy)(i,j) != 0.0)
				cout << "\nEnergy(" << i << ',' << j << ") = " << (*m_pairEnergy)(i,j);
}
#endif

} //namespace ogdf
