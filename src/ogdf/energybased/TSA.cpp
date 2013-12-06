/** \file
 * \brief  Iimplementation of class TSA
 *
 * This class realizes the TSA Algorithm for
 * automtatic graph drawing. It minimizes the energy
 * of the drawing using thermodynamic simulated annealing. This file
 * contains the main simulated annealing algorithm and
 * the fnction for computing the next candidate layout
 * that should be considered.
 *
 * \author Olivier Deckers
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

#include <ogdf/energybased/TSA.h>
#include <ogdf/basic/Math.h>
#include <time.h>

//TODO: in addition to the layout size, node sizes should be used in
//the initial radius computation in case of "all central" layouts with
//huge nodes
//the combinations for parameters should be checked: its only useful
//to have a slow shrinking if you have enough time to shrink down to
//small radius

namespace ogdf {

	const double TSA::m_startingTemp = 1e-2;
	const double TSA::m_defaultEndTemperature = 1e-5;

	//initializes internal data and the random number generator
	TSA::TSA():
	m_temperature(m_startingTemp),
	m_energy(0.0),
	m_endTemperature(m_defaultEndTemperature),
	m_quality(1.0)
	{
		
	}

	//allow resetting in between subsequent calls
	void TSA::initParameters()
	{
		m_diskRadius = computeDiskRadius(m_temperature);
		m_energy = 0.0;
		//m_numberOfIterations = 0; //is set in member function

		unsigned int t = (unsigned) time(NULL);
		//srand(t);
		cout << "seed: " << t << endl;
		srand(t);
		//srand((unsigned int) 1385114936);
	}

	void TSA::setQuality(double quality) 
	{
		OGDF_ASSERT(quality >= 0);
		m_quality = quality;
	}

	void TSA::setStartTemperature(double startTemp)
	{
		OGDF_ASSERT(startTemp >= 0);
		m_temperature=startTemp;
	}

	//whenever an energy function is added, the initial energy of the new function
	//is computed and added to the initial energy of the layout
	void TSA::addEnergyFunction(EnergyFunction *F, double weight)
	{
		m_energyFunctions.pushBack(F);
		OGDF_ASSERT(weight >= 0);
		m_weightsOfEnergyFunctions.pushBack(weight);
		F->computeEnergy();
		m_energy += F->energy();
	}

	List<string> TSA::returnEnergyFunctionNames()
	{
		List<string> names;
		ListIterator<EnergyFunction*> it;
		for(it = m_energyFunctions.begin(); it.valid(); it = it.succ())
			names.pushBack((*it)->getName());
		return names;
	}

	List<double> TSA::returnEnergyFunctionWeights()
	{
		List<double> weights;
		ListIterator<double> it;
		for(it = m_weightsOfEnergyFunctions.begin(); it.valid(); it = it.succ())
			weights.pushBack(*it);
		return weights;
	}

	//newVal is the energy value of a candidate layout. It is accepted if it is lower
	//than the previous energy of the layout or if m_fineTune is not tpFine and
	//the difference to the old energy divided by the temperature is smaller than a
	//random number between zero and one
	bool TSA::testEnergyValue(double newVal)
	{
		bool accepted = true;
		if(newVal > m_energy) {
			accepted = false;

			double testval = exp((m_energy-newVal)/ m_temperature);
			double compareVal = randNum(); // number between 0 and 1

			if(compareVal < testval)
				accepted = true;

		}
		return accepted;
	}

	//divides number returned by rand by RAND_MAX to get number between zero and one
	inline double TSA::randNum() const
	{
		double val = rand();
		val /= RAND_MAX;
		return val;
	}

	//chooses random vertex and a new random position for it on a circle with radius m_diskRadius
	//around its previous position
	edge TSA::computeCandidateLayout(
	const GraphAttributes &AG,
	DPoint &newSourcePos, DPoint &newTargetPos) const
	{
		int randomPos = randomNumber(0,m_nonIsolatedNodes.size()-1);
		node s = *(m_nonIsolatedNodes.get(randomPos));
		int randomNeighbour = randomNumber(0, s->degree()-1);
		adjEntry ae = s->firstAdj();
		for(int i=0; i<randomNeighbour; i++)
			ae = ae->succ();
		edge e = ae->theEdge();

		s = e->source();
		node t = e->target();

		double oldx = AG.x(s);
		double oldy = AG.y(s);
		double randomAngle = randNum() * 2.0 * Math::pi;
		newSourcePos.m_y = oldy+sin(randomAngle)*m_diskRadius * randNum();
		newSourcePos.m_x = oldx+cos(randomAngle)*m_diskRadius * randNum();

		oldx = AG.x(t);
		oldy = AG.y(t);
		randomAngle = randNum() * 2.0 * Math::pi;
		newTargetPos.m_y = oldy+sin(randomAngle)*m_diskRadius * randNum();
		newTargetPos.m_x = oldx+cos(randomAngle)*m_diskRadius * randNum();

		return e;
	}

	//steps through all energy functions and adds the initial energy computed by each
	//function for the start layout
	void TSA::computeInitialEnergy()
	{
		OGDF_ASSERT(!m_energyFunctions.empty());
		ListIterator<EnergyFunction*> it;
		ListIterator<double> it2;
		it2 = m_weightsOfEnergyFunctions.begin();
		for(it = m_energyFunctions.begin(); it.valid() && it2.valid(); it=it.succ(), it2 = it2.succ())
			m_energy += (*it)->energy() * (*it2);
	}

	//the vertices with degree zero are placed below all other vertices on a horizontal
	// line centered with repect to the rest of the drawing
	void TSA::placeIsolatedNodes(GraphAttributes &AG) const {
		double minX = 0.0;
		double minY = 0.0;
		double maxX = 0.0;
		double maxY = 0.0;

		if(!m_nonIsolatedNodes.empty()) {
			//compute a rectangle that includes all non-isolated vertices
			node v = m_nonIsolatedNodes.front();
			minX = AG.x(v);
			minY = AG.y(v);
			maxX = minX;
			maxY = minY;
			ListConstIterator<node> it;
			for(it = m_nonIsolatedNodes.begin(); it.valid(); ++it) {
				v = *it;
				double xVal = AG.x(v);
				double yVal = AG.y(v);
				double halfHeight = AG.height(v) / 2.0;
				double halfWidth = AG.width(v) / 2.0;
				if(xVal - halfWidth < minX) minX = xVal - halfWidth;
				if(xVal + halfWidth > maxX) maxX = xVal + halfWidth;
				if(yVal - halfHeight < minY) minY = yVal - halfHeight;
				if(yVal + halfHeight > maxY) maxY = yVal + halfHeight;
			}
		}

		// compute the width and height of the largest isolated node
		List<node> isolated;
		node v;
		const Graph &G = AG.constGraph();
		double maxWidth = 0;
		double maxHeight = 0;
		forall_nodes(v,G) if(v->degree() == 0) {
			isolated.pushBack(v);
			if(AG.height(v) > maxHeight) maxHeight = AG.height(v);
			if(AG.width(v) > maxWidth) maxWidth = AG.width(v);
		}
		// The nodes are placed on a line in the middle under the non isolated vertices.
		// Each node gets a box sized 2 maxWidth.
		double boxWidth = 2.0*maxWidth;
		double commonYCoord = minY-(1.5*maxHeight);
		double XCenterOfDrawing = minX + ((maxX-minX)/2.0);
		double startXCoord = XCenterOfDrawing - 0.5*(isolated.size()*boxWidth);
		ListIterator<node> it;
		double xcoord = startXCoord;
		for(it = isolated.begin(); it.valid(); ++it) {
			v = *it;
			AG.x(v) = xcoord;
			AG.y(v) = commonYCoord;
			xcoord += boxWidth;
		}
	}



	//this is the main optimization routine with the loop that lowers the temperature
	//and the disk radius geometrically until the temperature is zero. For each
	//temperature, a certain number of new positions for a random vertex are tried
	void TSA::call(GraphAttributes &AG)
	{
		initParameters();

		OGDF_ASSERT(!m_energyFunctions.empty());

		const Graph &G = AG.constGraph();
		//compute the list of vertices with degree greater than zero
		G.allNodes(m_nonIsolatedNodes);
		ListIterator<node> it,itSucc;
		for(it = m_nonIsolatedNodes.begin(); it.valid(); it = itSucc) {
			itSucc = it.succ();
			if((*it)->degree() == 0) m_nonIsolatedNodes.del(it);
		}


		if(G.numberOfEdges() > 0) { //else only isolated nodes
			computeInitialEnergy();
			
			double totalCostDiff = 0, totalEntropyDiff = 0;
			double costDiff;
			int i = 0;
			int iterationsSinceLastChange = 0;
			//this is the main optimization loop
			while((m_temperature > m_endTemperature || i < 20) && m_diskRadius >= 1) {

				DPoint newSourcePos, newTargetPos;
				//choose random vertex and new position for vertex
				edge e = computeCandidateLayout(AG,newSourcePos,newTargetPos);

				//compute candidate energy and decide if new layout is chosen
				ListIterator<EnergyFunction*> it;
				ListIterator<double> it2 = m_weightsOfEnergyFunctions.begin();
				double newEnergy = 0.0;
				for(it = m_energyFunctions.begin(); it.valid(); it = it.succ()) {
					newEnergy += (*it)->computeCandidateEnergy(e,newSourcePos,newTargetPos) * (*it2);
					it2 = it2.succ();
				}
				OGDF_ASSERT(newEnergy >= 0.0);

				costDiff = newEnergy - m_energy;

				//this tests if the new layout is accepted. If this is the case,
				//all energy functions are informed that the new layout is accepted
				if(testEnergyValue(newEnergy)) {
					totalCostDiff += costDiff;

					for(it = m_energyFunctions.begin(); it.valid(); it = it.succ())
						(*it)->candidateTaken();

					AG.x(e->source()) = newSourcePos.m_x;
					AG.y(e->source()) = newSourcePos.m_y;
					AG.x(e->target()) = newTargetPos.m_x;
					AG.y(e->target()) = newTargetPos.m_y;
					m_energy = newEnergy;

					iterationsSinceLastChange = 0;
				}

				if(costDiff > 0) {
					totalEntropyDiff -= costDiff / m_temperature;
				}
				
				if(totalCostDiff >= 0 || totalEntropyDiff == 0) {
					m_temperature = m_startingTemp;
				}
				else {
					m_temperature = m_quality * (totalCostDiff / totalEntropyDiff);
				}

				m_diskRadius = computeDiskRadius(m_temperature);

				cout << "temperature: " << m_temperature << endl;
				cout << "diskradius: " << m_diskRadius << endl;
				cout << "energy: " << m_energy << endl;
				cout << "iteration: " << i << endl;
				iterationsSinceLastChange ++;

				if(iterationsSinceLastChange > 1000)
					break;

				i ++;
			}
		}
		//if there are zero degree vertices, they are placed using placeIsolatedNodes
		if(m_nonIsolatedNodes.size() != G.numberOfNodes())
			placeIsolatedNodes(AG);

	}

	double TSA::computeDiskRadius(double temperature) const {
		return min(500.0, 1e4 * (temperature - 1e-6)); //TODO enhance
	}
} //namespace
