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
import java.util.Set;

/**
 *
 * @author Alfonso E. Romero
 */
public abstract class GraphSimilarity {

    Assignment annotations;

    public GraphSimilarity(Assignment annotations) {
        this.annotations = annotations;
    }

    public abstract float similarity(Set<GOTerm> s1, Set<GOTerm> s2);

    public abstract void setMaxAnnotations(double annot);
}
