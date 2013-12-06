/*
 * $Revision: 2565 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2012-07-07 17:14:54 +0200 (Sa, 07. Jul 2012) $
 ***************************************************************/

/** \file
 * \brief Implementation of class Planarity
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

#include <ogdf/internal/energybased/Planarity.h>

namespace ogdf {

	Planarity::~Planarity()
	{
		delete m_edgeNums;
		delete m_crossingMatrix;
	}


	// intializes number of edges and allocates memory for  crossingMatrix
	Planarity::Planarity(GraphAttributes &AG):
	EnergyFunction("Planarity",AG)
	{
		m_edgeNums = OGDF_NEW EdgeArray<int>(m_G,0);
		m_G.allEdges(m_nonSelfLoops);
		ListIterator<edge> it, itSucc;
		for(it = m_nonSelfLoops.begin(); it.valid(); it = itSucc) {
			itSucc = it.succ();
			if((*it)->isSelfLoop()) m_nonSelfLoops.del(it);
		}
		int e_num = 1;
		for(it = m_nonSelfLoops.begin(); it.valid(); ++it) (*m_edgeNums)[*it] = e_num ++;
		e_num --;
		m_crossingMatrix = new Array2D<double> (1,e_num,1,e_num);
	}


	// computes energy of layout, stores it and sets the crossingMatrix
	void Planarity::computeEnergy()
	{
		int e_num = m_nonSelfLoops.size();
		double energySum = 0;
		Array<edge> numEdge(1,e_num);
		edge e;
		ListIterator<edge> it;

		for(it = m_nonSelfLoops.begin(); it.valid(); ++it)
			numEdge[(*m_edgeNums)[*it]] = *it;
		for(int i = 1; i < e_num; i++) {
			e = numEdge[i];
			for(int j = i+1; j <= e_num ; j++) {
				double energy = 0;
				bool cross = intersect(e,numEdge[j], energy);
				(*m_crossingMatrix)(i,j) = cross ? energy : 0;
				if(cross) energySum += energy;
			}
		}
		m_energy = energySum;
	}


	// tests if two edges cross
	bool Planarity::intersect(const edge e1, const edge e2, double &energy) const
	{
		node v1s = e1->source();
		node v1t = e1->target();
		node v2s = e2->source();
		node v2t = e2->target();

		bool cross = false;
		DPoint inter;
		if(v1s != v2s && v1s != v2t && v1t != v2s && v1t != v2t)
			cross = lowLevelIntersect(currentPos(v1s),currentPos(v1t), currentPos(v2s),currentPos(v2t), energy);
		return cross;
	}


	// tests if two lines given by four points cross
	bool Planarity::lowLevelIntersect(
		const DPoint &e1s,
		const DPoint &e1t,
		const DPoint &e2s,
		const DPoint &e2t,
		double &energy) const
	{
		DPoint s1(e1s),t1(e1t),s2(e2s),t2(e2t);
		DLine l1(s1,t1), l2(s2,t2);
		DPoint dummy;
		bool intersect = l1.intersection(l2,dummy);
		if(intersect) {
			double length = l1.length()/2.0;
			double interDist = min(t1.distance(dummy), s1.distance(dummy));

			energy = interDist / length;
		}
		else {
			energy = 0;
		}
		return intersect;
	}


	void Planarity::compCandEnergy()
	{
		m_crossingChanges.clear();
		m_candidateEnergy = energy();

		node s = sourceTestNode();
		node t = targetTestNode();

		compCandEnergy(s, t, sourceTestPos());
		
		if(t != NULL) {
			compCandEnergy(t, s, targetTestPos());

			//calculate candidate energy for moved edge
			edge e,f;
			forall_adj_edges(e, s) if(e->target() == t) {
				int e_num = (*m_edgeNums)[e];
				ListIterator<edge> it;
				for(it = m_nonSelfLoops.begin(); it.valid(); ++it) if(*it != e) {
					f = *it;
					node s2 = f->source();
					node t2 = f->target();
					int f_num = (*m_edgeNums)[f];

					double intersectEnergy = 0;
					bool cross = lowLevelIntersect(sourceTestPos(),targetTestPos(),currentPos(s2),currentPos(t2), intersectEnergy);
					double priorIntersectEnergy = (*m_crossingMatrix)(min(e_num,f_num),max(e_num,f_num));

					if(priorIntersectEnergy != intersectEnergy) {
						m_candidateEnergy -= priorIntersectEnergy; // this intersection was saved
						m_candidateEnergy += intersectEnergy; // produced a new intersection
						ChangedCrossing cc;
						cc.edgeNum1 = min(e_num,f_num);
						cc.edgeNum2 = max(e_num,f_num);
						cc.crossEnergy = intersectEnergy;
						m_crossingChanges.pushBack(cc);
					}
				}
			}
		}
	}

	// computes the energy if the node returned by testNode() is moved
	// to position testPos().
	void Planarity::compCandEnergy(const node v, const node ignore, const DPoint newPos)
	{
		edge e;

		forall_adj_edges(e,v) if(!e->isSelfLoop()) {
			// first we compute the two endpoints of e if v is on its new position
			node s = e->source();
			node t = e->target();

			if((s == v && t == ignore) || (s == ignore && t == v)) {
				continue;
			}

			DPoint p1 = newPos;
			DPoint p2 = (s==v)? currentPos(t) : currentPos(s);
			int e_num = (*m_edgeNums)[e];
			edge f;
			// now we compute the crossings of all other edges with e
			ListIterator<edge> it;
			for(it = m_nonSelfLoops.begin(); it.valid(); ++it) if(*it != e) {
				f = *it;
				node s2 = f->source();
				node t2 = f->target();

				if((s2 == v && t2 == ignore) || (s2 == ignore && t2 == v)) {
					continue;
				}

				if(s2 != s && s2 != t && t2 != s && t2 != t) {
					double intersectEnergy = 0;
					bool cross = lowLevelIntersect(p1,p2,currentPos(s2),currentPos(t2), intersectEnergy);
					int f_num = (*m_edgeNums)[f];
					double priorIntersectEnergy = (*m_crossingMatrix)(min(e_num,f_num),max(e_num,f_num));
					
					if(priorIntersectEnergy != intersectEnergy) {
						m_candidateEnergy -= priorIntersectEnergy; // this intersection was saved
						m_candidateEnergy += intersectEnergy; // produced a new intersection
						ChangedCrossing cc;
						cc.edgeNum1 = min(e_num,f_num);
						cc.edgeNum2 = max(e_num,f_num);
						cc.crossEnergy = intersectEnergy;
						m_crossingChanges.pushBack(cc);
					}
				}
			}
		}
	}


	// this functions sets the crossingMatrix according to candidateCrossings
	void Planarity::internalCandidateTaken() {
		ListIterator<ChangedCrossing> it;
		for(it = m_crossingChanges.begin(); it.valid(); ++ it) {
			ChangedCrossing cc = *(it);
			(*m_crossingMatrix)(cc.edgeNum1,cc.edgeNum2) = cc.crossEnergy;
		}
	}


#ifdef OGDF_DEBUG
void Planarity::printInternalData() const {
	cout << "\nCrossing Matrix:";
	int e_num = m_nonSelfLoops.size();
	for(int i = 1; i < e_num; i++) {
		cout << "\n Edge " << i << " crosses: ";
		for(int j = i+1; j <= e_num; j++)
			if((*m_crossingMatrix)(i,j)) cout << j << " ";
	}
	cout << "\nChanged crossings:";
	if(sourceTestNode() == NULL) cout << " None.";
	else {
		ListConstIterator<ChangedCrossing> it;
		for(it = m_crossingChanges.begin(); it.valid(); ++it) {
			ChangedCrossing cc = *(it);
			cout << " (" << cc.edgeNum1 << "," << cc.edgeNum2 << ")" << cc.crossEnergy;
		}
	}
}
#endif

}
