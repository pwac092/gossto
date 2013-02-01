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

import java.util.*;

/**
 *
 * @author aeromerow
 */
public class GOTerm implements Comparable<GOTerm> {

    /**
     * Identifier of the GO term, for example GO:0000187
     */
    private String id;
    /**
     * The integer identifier to the term, in the previous case 187
     */
    private int numId;
    /**
     * Description of the function
     */
    private String function;
    /**
     * Cache of ancestors
     */
    private Set<GOTerm> ancestorsCache;
    /**
     * More general terms this belong to, indexed by identifier of the relations
     */
    private Map<String, List<GOTerm>> parents;
    /**
     * More specific terms belonging to this, indexed by identifier of the
     * relations
     */
    private Map<String, List<GOTerm>> children;
    /**
     * Set of relations indexed
     */
    private static Set<String> relations;
    /**
     * Sub-ontology it belongs to
     */
    Ontology ontology;

    public GOTerm(final String GOid, final String function) throws IllegalArgumentException {
        this.id = GOid;
        this.function = function;
        this.parents = new HashMap<String, List<GOTerm>>();
        this.children = new HashMap<String, List<GOTerm>>();
        this.ontology = null;

        try {
            this.numId = Integer.parseInt(GOid.split(":")[1]);
        } catch (ArrayIndexOutOfBoundsException ex) {
            throw new IllegalArgumentException("GO term identifier " + GOid + " is not valid");
        }

        if (relations == null) {
            throw new NullPointerException("ERROR: The kind of relations that are going to be parsed should be set first via GOTerm.setRelations");
        }

        for (String rel : relations) {
            this.parents.put(rel, new ArrayList<GOTerm>());
            this.children.put(rel, new ArrayList<GOTerm>());
        }
        this.ancestorsCache = null;
    }

    public GOTerm(String GOid, String function, Ontology onto) throws IllegalArgumentException {
        this(GOid, function);
        this.ontology = onto;
        this.ancestorsCache = null;
    }

    public void setOntology(Ontology ontology) {
        this.ontology = ontology;
    }

    /**
     * Return the set of ancestors (including itself) for a certain given
     * relation
     *
     * @param rel string with the relation
     * @return set of GOTerms ancestors of this
     */
    public Set<GOTerm> getAncestorsRelation(final String rel) {
        Set<GOTerm> ret = new HashSet<GOTerm>();
        Queue<GOTerm> queue = new LinkedList<GOTerm>();
        queue.add(this);

        while (!queue.isEmpty()) {
            final GOTerm term = queue.poll();

            if (!ret.contains(term)) {
                ret.add(term);
                queue.addAll(term.getParentsForRelation(rel));
            }
        }
        return ret;
    }

    /**
     * Return the set of all ancestors (including itself) for all indexed
     * relations
     *
     * @return set of GOTerms ancestors of this
     */
    public Set<GOTerm> getAncestors() {

        if (this.ancestorsCache == null) {
            Set<GOTerm> ret = new HashSet<GOTerm>();
            Queue<GOTerm> queue = new LinkedList<GOTerm>();
            queue.add(this);

            while (!queue.isEmpty()) {
                final GOTerm term = queue.poll();

                if (!ret.contains(term)
                        && term.ontology == this.ontology) {
                    ret.add(term);
                    for (String rel : GOTerm.relations) {
                        queue.addAll(term.getParentsForRelation(rel));
                    }
                }
            }
            this.ancestorsCache = ret;
            return ret;
        } else {
            return this.ancestorsCache;
        }
    }

    /**
     * Return the set of all descendants (including itself) for all indexed
     * relations
     *
     * @return set of GOTerms descendants of this
     */
    public Set<GOTerm> getDescendants() {

        Set<GOTerm> ret = new HashSet<GOTerm>();
        Queue<GOTerm> queue = new LinkedList<GOTerm>();
        queue.add(this);

        while (!queue.isEmpty()) {
            final GOTerm term = queue.poll();

            if (!ret.contains(term)) {
                ret.add(term);
                for (String rel : GOTerm.relations) {
                    queue.addAll(term.getChildrenForRelation(rel));
                }
            }
        }
        return ret;
    }

    public void addParent(final String relation, final GOTerm term) {
        if (checkRelation(relation)) {
            // only parents from the same ontology are considered
            if (term.ontology == this.ontology) {
                this.parents.get(relation).add(term);
            }
        } else {
            throw new IllegalArgumentException("ERROR: The relation " + relation + " is not being"
                    + " parsed in this instance of the Gene Ontology");
        }
    }

    public void addChildren(final String relation, final GOTerm term) {
        if (checkRelation(relation)) {
            // only children from the same ontology are considered
            if (term.ontology == this.ontology) {
                this.children.get(relation).add(term);
            }
        } else {
            throw new IllegalArgumentException("ERROR: The relation " + relation + " is not being"
                    + " parsed in this instance of the Gene Ontology");
        }
    }

    /**
     * Return the set of children (including itself) for a certain given
     * relation
     *
     * @param rel string with the relation
     * @return set of GOTerms children of this
     */
    public Set<GOTerm> getDescendantsRelation(final String rel) {
        Set<GOTerm> ret = new HashSet<GOTerm>();
        Queue<GOTerm> queue = new LinkedList<GOTerm>();
        queue.add(this);

        while (!queue.isEmpty()) {
            GOTerm term = queue.poll();

            if (!ret.contains(term)) {
                ret.add(term);
                List<GOTerm> yeah = term.getChildrenForRelation(rel);
                queue.addAll(yeah);
            }
        }
        return ret;
    }

    public static void setRelations(String[] rel) {
        relations = new HashSet<String>(Arrays.asList(rel));
    }

    public boolean isRoot(final String rel) {
        return this.getParentsForRelation(rel).isEmpty();
    }

    public boolean isLeaf(final String rel) {
        return this.getChildrenForRelation(rel).isEmpty();
    }

 
    public Ontology getOntology() {
        return ontology;
    }

    public int getNumericId() {
        return this.numId;
    }

    public String getGOid() {
        return this.id;
    }

    public void setGOid(final String GOid) {
        this.id = GOid;
    }

    public String getFunction() {
        return function;
    }

    public void setFunction(String function) {
        this.function = function;
    }

    public final List<GOTerm> getChildrenForRelation(final String rel) {
        if (checkRelation(rel)) {
            return this.children.get(rel);
        } else {
            throw new IllegalArgumentException("ERROR: The relation " + rel + " has not been parsed in this instance of the Gene Ontology");
        }
    }

    public final List<GOTerm> getParentsForRelation(final String rel) {
        if (checkRelation(rel)) {
            return this.parents.get(rel);
        } else {
            throw new IllegalArgumentException("ERROR: The relation " + rel + " has not been parsed in this instance of the Gene Ontology");
        }
    }

    public static boolean checkRelation(final String rel) {
        return relations != null ? relations.contains(rel) : false;
        /*if (relations != null) {

         return ;
         } else {
         return false;
         }*/
    }

    @Override
    public int hashCode() {
        return this.numId;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final GOTerm other = (GOTerm) obj;
        if (this.numId != other.numId) {
            return false;
        }
        return true;
    }

    @Override
    public int compareTo(GOTerm t) {
        return this.getNumericId() - t.getNumericId();
    }
}
