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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 *
 * @author aeromero
 */
public class GeneOntology { //extends Singleton {

    /** Makes a correspondence between every ontology name and itself */
    Map<String, Ontology> ontologies;
    /** Maps strings to (non obsolete) GO terms */
    Map<String, GOTerm> goIdToTerm;
    /** Map identifiers of obsolete GO terms to their alternatives (marked as
    "consider_using:" */
    Map<String, List<GOTerm>> dummyTermsToGOTerms;

    public GeneOntology() {
        this.ontologies = new HashMap<String, Ontology>();
        this.goIdToTerm = new HashMap<String, GOTerm>();
        this.dummyTermsToGOTerms = new HashMap<String, List<GOTerm>>();
    }

    public boolean isNonObsolete(String id) {
        return this.goIdToTerm.containsKey(id);
    }

    public GOTerm getTermById(String id) {
        return this.goIdToTerm.get(id);
    }

    public void putTermById(GOTerm term, String id) {
        this.goIdToTerm.put(id, term);
    }

    public List<GOTerm> getSynonyms(String goTermId) {
        if (this.dummyTermsToGOTerms.containsKey(goTermId)) {
            return this.dummyTermsToGOTerms.get(goTermId);
        } else {
            return new ArrayList<GOTerm>();
        }
    }

    public void addSynonymForTerm(String obsoleteTerm, String alternativeTerm) {
        if (this.isNonObsolete(alternativeTerm)) {

            GOTerm alternative = this.goIdToTerm.get(alternativeTerm);

            if (this.dummyTermsToGOTerms.containsKey(obsoleteTerm)) {
                List<GOTerm> alternatives = this.dummyTermsToGOTerms.get(obsoleteTerm);
                if (!alternatives.contains(alternative)) {
                    alternatives.add(alternative);
                }
            } else {
                List<GOTerm> alternatives = new ArrayList<GOTerm>();
                alternatives.add(alternative);
                this.dummyTermsToGOTerms.put(obsoleteTerm, alternatives);
            }
        } else {

            // The term has not been found, we test several levels of indirection
            if (this.dummyTermsToGOTerms.containsKey(alternativeTerm)) {
                List<GOTerm> alternativesToTheAlternative = this.dummyTermsToGOTerms.get(alternativeTerm);
                this.dummyTermsToGOTerms.put(obsoleteTerm, alternativesToTheAlternative);
            } else {
                // this is a strange case: it means that a term is obsolete, and points to another
                // term which is obsolete and either still no alternative terms have been
                // processed for it, or, there is more than a level of indirection
                throw new RuntimeException("ERROR: Too many indirection levels found for GO term id: " + obsoleteTerm);
            }
        }
    }

    /*
    public void addSynonymForTerm(String originalId, String synonymId) {
    if (this.existsId(originalId)) {
    GOTerm term = this.getTermById(originalId);

    if (!this.existsId(synonymId)) {
    this.putTermById(term, synonymId);
    }
    }
    }
     */
    public boolean hasOntology(String ontologyName) {
        return this.ontologies.containsKey(ontologyName);
    }

    public Ontology getOntology(String ontologyName) {
        return this.ontologies.get(ontologyName);
    }

    public void addTermToOntology(GOTerm term, String ontologyName) {
        if (this.hasOntology(ontologyName)) {
            Ontology onto = this.getOntology(ontologyName);
            onto.addNode(term);
            term.setOntology(onto);
        } else {
            Ontology onto = new Ontology(ontologyName);
            onto.addNode(term);
            this.ontologies.put(ontologyName, onto);
            term.setOntology(onto);
        }
    }

    public void putTermByListOfIds(GOTerm term, List<String> l) {
        for (String id : l) {
            this.putTermById(term, id);
        }
    }
}
