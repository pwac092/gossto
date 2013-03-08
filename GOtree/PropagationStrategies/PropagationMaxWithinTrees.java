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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package GOtree.PropagationStrategies;

import GOtree.GOTerm;
import GOtree.GeneOntology;
import GOtree.PropagationStrategy; 
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;
import java.util.Queue;

/**
 *
 * @author aeromero
 */
public class PropagationMaxWithinTrees extends PropagationStrategy {

    @Override
    public Map<String, Double> propagateRow(final Map<String, Double> inputValues, final GeneOntology ontology) {

        Map<String, Double> outputValues = new HashMap<String, Double>();

        for (Map.Entry<String, Double> inputPair : inputValues.entrySet()) {

            String termIdentifier = inputPair.getKey();
            double value = inputPair.getValue();

            Queue<GOTerm> terms = new LinkedList<GOTerm>();

            if (ontology.isNonObsolete(termIdentifier)) {
                GOTerm _term = ontology.getTermById(termIdentifier);
                terms.add(_term);
            } else {
                terms.addAll(ontology.getSynonyms(termIdentifier));
            }

            while (!terms.isEmpty()) {

                GOTerm term = terms.poll();

                for (GOTerm propagatedTerm : term.getAncestors()) {

                    if (propagatedTerm.getOntology() == term.getOntology()) { // propagate only within the same ontology

                        String goTermId = propagatedTerm.getGOid();
                        if (!outputValues.containsKey(goTermId)) {
                            outputValues.put(goTermId, value);
                        } else {
                            double oldValue = outputValues.get(goTermId);
                            outputValues.put(goTermId, Math.max(oldValue, value));
                        }
                    }
                }
            }
        }

        return outputValues;
    }
}
