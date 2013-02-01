/*This file is part of GOSSTO.
    GOSSTO is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GOSSTO is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GOSSTO.  If not, see <http://www.gnu.org/licenses/>.
*/
package GOtree;



import java.io.BufferedReader; 
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author aeromero
 */
public class GeneOntologyParser {

    Map<String, List<String>> synonym;
    Map<String, List<String>> is_a;
    Map<String, Map<String, List<String>>> relations;
    String goId = "", name = "", ontology = "";
    List<String> parents = new ArrayList<String>();
    List<String> consider = new ArrayList<String>();
    List<String> alternative_ids = new ArrayList<String>();
    Map<String, List<String>> related = new HashMap<String, List<String>>();
    GeneOntology result;

    public GeneOntologyParser() {
        synonym = new HashMap<String, List<String>>();
        is_a = new HashMap<String, List<String>>();
        relations = new HashMap<String, Map<String, List<String>>>();
        result = new GeneOntology();
    }

    private void clear() {
        this.synonym.clear();
        this.is_a.clear();
        this.relations.clear();
    }

    private void clearForTerm() {
        parents = new ArrayList<String>();
        consider = new ArrayList<String>();
        related = new HashMap<String, List<String>>();
        alternative_ids = new ArrayList<String>();
        goId = "";
        name = "";
        ontology = "";
    }
    
    public GeneOntology readFromOBOFile(String OBOfileName) throws IOException, FileNotFoundException
    {
        return this.readFromOBOFile(OBOfileName, true);
    }

    public GeneOntology readFromOBOFile(String OBOfileName, boolean useConsiderAsSynonym) throws IOException, FileNotFoundException {
        this.clear();

        BufferedReader reader = new BufferedReader(new FileReader(OBOfileName));

        parents = new ArrayList<String>();
        consider = new ArrayList<String>();
        related = new HashMap<String, List<String>>();

        boolean isObsolete = false, processingATerm = false;
        boolean stillInHeader = true, locked = false;

        while (reader.ready()) {
            String line = reader.readLine().trim();
           
            // we discard non useful information
            if (line.startsWith("xref") || line.startsWith("subset") || line.startsWith("synonym")
                    || line.startsWith("created_by:") || line.startsWith("creation_date")
                    || line.startsWith("comment:") || line.startsWith("intersection_of")
                    || line.startsWith("def:")) {
                continue;
            }

            if (line.startsWith("[Term]")) {
                if (!stillInHeader) {
                    // process term and store it
                    if (isObsolete) {
                        if (useConsiderAsSynonym && !consider.isEmpty()) {
                            this.synonym.put(goId, consider);
                        }
                        isObsolete = false;
                    } else if (processingATerm) {
                        processGOterm(goId, name, ontology);
                    }
                    this.clearForTerm();
                }
                processingATerm = true;
                stillInHeader = false;

            } else if (line.startsWith("[Typedef]")) {
                if (isObsolete) {
                    if (useConsiderAsSynonym && !consider.isEmpty()) {
                        this.synonym.put(goId, consider);
                    }
                    isObsolete = false;
                } else if (processingATerm) {
                    processGOterm(goId, name, ontology);
                }

                locked = true;
                isObsolete = false;
                processingATerm = false;
                this.clearForTerm();

            } else if (!stillInHeader && !locked) {

                if (line.startsWith("id:")) {
                    goId = (line.split("\\s+")[1]).trim();
                } else if (line.startsWith("name:")) {
                    name = line.split("\\s+", 2)[1].trim();
                } else if (line.startsWith("alt_id:")) {
                    alternative_ids.add(line.split("\\s+")[1].trim());
                } else if (line.startsWith("namespace:")) {
                    ontology = (line.split("\\s+")[1]).trim();
                } else if (line.startsWith("is_a:")) {
                    parents.add((line.split("\\s+")[1]).trim());
                } else if (line.startsWith("is_obsolete:")) {
                    isObsolete = true;
                } else if (line.startsWith("consider:")) {
                    consider.add(line.split("\\s+")[1].trim());
                } else if (line.startsWith("relationship:")) {
                    String fields[] = line.split("\\s+");
                    String type = fields[1].trim();
                    String relatedGOTerm = fields[2].trim();
                    if (!this.related.containsKey(type)) {
                        List<String> values = new LinkedList<String>();
                        values.add(relatedGOTerm);
                        this.related.put(type, values);
                    } else {
                        this.related.get(type).add(relatedGOTerm);
                    }
                }
            }

            locked = (locked && line.length() == 0);
        }

        if (!isObsolete && processingATerm) {
            processGOterm(goId, name, ontology);
        }

        this.adjustRelations();
        reader.close();
        return result;
    }

    private void processGOterm(String goId, String name, String ontology) {
        GOTerm term = new GOTerm(goId, name);
        result.addTermToOntology(term, ontology);
        result.putTermById(term, goId);
        is_a.put(goId, parents);

        if (!this.alternative_ids.isEmpty())
        {
            for (String alternative_id : alternative_ids)
                result.putTermById(term, alternative_id);
        }

        if (!this.related.isEmpty()) {
            this.relations.put(goId, related);
        }
        if (!consider.isEmpty()) {
            System.err.println("ERROR: consider tag present for: " + goId);
            this.synonym.put(goId, consider);
        }
    }

    private void adjustRelations() {

        for (Map.Entry<String, List<String>> pair : this.is_a.entrySet()) {
            GOTerm term = result.getTermById(pair.getKey());
            for (String parent : pair.getValue()) {
                GOTerm parent_term = result.getTermById(parent);
                term.addParent("is_a", parent_term);
                parent_term.addChildren("is_a", term);
            }
        }

        for (String baseTerm : this.relations.keySet()) {

            for (String relation : this.relations.get(baseTerm).keySet()) {

                if (GOTerm.checkRelation(relation)) {
                    GOTerm base = result.getTermById(baseTerm);

                    for (String relatedTerm : this.relations.get(baseTerm).get(relation)) {
                        GOTerm relatedGOTerm = result.getTermById(relatedTerm);
                        base.addParent(relation, relatedGOTerm);
                        relatedGOTerm.addChildren(relation, base);
                    }
                }
            }
        }

        for (String obsoleteTerm : synonym.keySet()) {
            for (String alternativeTerm : synonym.get(obsoleteTerm)) {
                result.addSynonymForTerm(obsoleteTerm, alternativeTerm);
            }
        }
    }
}
