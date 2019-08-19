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

 


import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 *
 * @author aeromero
 */
public class Ontology {

    /** name of the ontology */
    String name;
    /** Set of GO terms of this ontology (no synonyms are allowed here), just
    effective GO terms */
    Map<String, GOTerm> goTermByID;
    /** Number of effective GO terms */
    int numGOterms;

    /** Main constructor
     *  @param name name of the ontology
     */
    public Ontology(String name) {
        this.numGOterms = 0;
        this.goTermByID = new HashMap<String, GOTerm>();
        this.name = name;
    }

    /** Returns true if the go identifier is in the ontology 
     *  @param GOid string with the GO identifier we want to check
     *  @return true if the identifier is in the ontology 
     */
    public boolean containsGOid(String GOid) {
        return this.goTermByID.containsKey(GOid);
    }

    /** Returns the GO term object corresponding to the identifier passed
     *  @param id identifier which is passed
     *  @return the object representing that GO Term
     */
    public GOTerm getTermByID(String id) {
        return this.goTermByID.get(id);
    }

    /** Returns the name of the ontology 
     *  @return name of the ontology 
     */
    public String getName() {
        return name;
    }

    /** Return the size of the ontology (number of terms)
     */
    public int size() {
        return this.goTermByID.size();
    }

    /** Adds a new node in the ontology
     *  @param node node we want to add in the ontology
     */
    public void addNode(GOTerm node) {
        this.goTermByID.put(node.getGOid(), node);
        this.numGOterms = this.goTermByID.size();
    }

    public List<GOTerm> getSetOfTerms() {
        return new ArrayList<GOTerm>(this.goTermByID.values());
    }
}
