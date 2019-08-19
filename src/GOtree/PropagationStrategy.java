/*This file is part of GOssTo.
    GOssTo is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GOssTo is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GOssTo.  If not, see <http://www.gnu.org/licenses/>.
*/

package GOtree;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor. 
 */


import java.util.Map;

/**
 *
 * @author aeromero
 */
abstract public class PropagationStrategy {

    protected abstract Map<String, Double> propagateRow(Map<String, Double> values, GeneOntology ontology);

    public Assignment propagateAssignment(Assignment other, GeneOntology onto) {
        Assignment propagated = new Assignment();
        for (String protein : other.getRowIdentifiers()) {
            Map<String, Double> propagatedRow = propagateRow(other.getGOTermScoresForProteinId(protein), onto);

            for (Map.Entry<String, Double> pairs : propagatedRow.entrySet()) {
                propagated.setValue(protein, pairs.getKey(), pairs.getValue());
            }
        }
        return propagated;
    }
}
