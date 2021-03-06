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
package HSM.GraphSimilarities;

import GOtree.Assignment;
import GOtree.GOTerm;
import java.util.HashSet;
import java.util.Set;

/**
 *
 * @author Alfonso E. Romero
 */
public class SimUISimilarity extends GraphSimilarity {

    public SimUISimilarity(Assignment annotations) {
        super(annotations);
    }
    
    @Override
    public float similarity(Set<GOTerm> s1, Set<GOTerm> s2) {
        Set<GOTerm> union = new HashSet<GOTerm>(s1);
        union.addAll(s2);
        Set<GOTerm> intersection = new HashSet<GOTerm>(s1);
        intersection.retainAll(s2);
        if (!union.isEmpty()) {
            return (float) intersection.size() / (float) union.size();
        } else {
            return 0.0f;
        }
    }

    @Override
    public void setMaxAnnotations(double annot) {
        //
    }
    
}
