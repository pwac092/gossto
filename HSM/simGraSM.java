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
package HSM;

import GOtree.Assignment;
import GOtree.GOTerm;
import Jama.Matrix;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import util.TinyLogger;

class TermWithIC implements Comparable<TermWithIC> {

    public GOTerm term;
    public float IC;

    public TermWithIC(GOTerm term, float IC) {
        this.term = term;
        this.IC = IC;
    }

    @Override
    public int compareTo(TermWithIC t) {
        if (this.IC < t.IC) {
            return -1;
        } else if (this.IC > t.IC) {
            return 1;
        } else {
            return 0;
        }
    }
}

/**
 *
 * @author Alfonso E. Romero
 */
//Implements the simGraSM semantic similarity measure
public class simGraSM extends HSM {

    /**
     * information content cache
     */
    Map<GOTerm, Float> icCache;

    /**
     * Caches the number of paths between two GO terms, using their goNumericId
     * (and forcing the first one to be less or equal than the second)
     */
    Map<Integer, Map<Integer, Integer>> numPathsCache;

    public simGraSM(GOTerm[] allTerms, String[] genes, String[][] goIds, GOTerm[][] axis, String[] targets, Assignment annotations, String[] relations, TinyLogger logw) {
        super(allTerms, genes, goIds, axis, targets, annotations, relations, logw);
        icCache = new HashMap<GOTerm, Float>();
        numPathsCache = new HashMap<Integer, Map<Integer, Integer>>();
    }

    //Returns the common ancestors of the GO terms 'target1' & 'target2'
    private HashSet<GOTerm> commonAncestors(GOTerm target1, GOTerm target2) {
        HashSet<GOTerm> CA = new HashSet<GOTerm>(target1.getAncestors()); //Common Ancestors
        CA.retainAll(target2.getAncestors());
        return CA;
    }

    //Recursive method to find paths from a given GO term 'target' to another 'ancestor'. 'path' is the current path, 'allPaths' are all found paths between the terms
    private void pathFinder(GOTerm ancestor, GOTerm currentNode, Set<GOTerm> path, Set<Set<GOTerm>> allPaths) {
        for (String rel : this.relations) {
            Set<GOTerm> parents = new HashSet<GOTerm>(currentNode.getParentsForRelation(rel));
            if (parents.contains(ancestor)) //If a path is found from target to ancestor, add it to the path list
            {
                allPaths.add(path);
                return;
            } else if (currentNode.isRoot(rel)) //If the current node is the root of an ontology, a path has not been found
            {
                return;
            }
            for (GOTerm parent : parents) {
                Set<GOTerm> myPath = new HashSet<GOTerm>(path);
                myPath.add(parent);
                pathFinder(ancestor, parent, myPath, allPaths);
            }
        }
    }

    private int getNumPathsTermAncestor(GOTerm term, GOTerm ancestor) {
        Set<GOTerm> path = new HashSet<GOTerm>();
        Set<Set<GOTerm>> allPaths = new HashSet<Set<GOTerm>>();
        pathFinder(ancestor, term, path, allPaths);
        return allPaths.size();
    }

    private int getNumPaths(GOTerm _t1, GOTerm _t2) {
        int t1 = Math.min(_t1.getNumericId(), _t2.getNumericId());
        int t2 = Math.max(_t1.getNumericId(), _t2.getNumericId());

        if (!numPathsCache.containsKey(t1) || !numPathsCache.get(t1).containsKey(t2)) {
            int numPaths;

            if (_t1.getAncestors().contains(_t2)) {
                // _t2 is an ancestor of _t1 (base)
                numPaths = getNumPathsTermAncestor(_t1, _t2);

            } else if (_t2.getAncestors().contains(_t1)) {
                // _t1 is an ancestor of _t2 (base)
                numPaths = getNumPathsTermAncestor(_t2, _t1);
            } else {
                // there is no path from t1 to t2 or vice versa
                numPaths = 0;
            }

            if (!numPathsCache.containsKey(t1)) {
                numPathsCache.put(t1, new HashMap<Integer, Integer>());
            }
            numPathsCache.get(t1).put(t2, numPaths);
            return numPaths;
        } else {
            return numPathsCache.get(t1).get(t2);
        }
    }

    private boolean DisjAnc(GOTerm c, GOTerm a1, GOTerm a2) {
        int nPaths = getNumPaths(a1, a2);
        int nPaths1 = getNumPaths(a1, c);
        int nPaths2 = getNumPaths(a2, c);
        return nPaths1 >= nPaths + nPaths2;
    }

    private Map<GOTerm, Float> computeIC(Set<GOTerm> CA, double numAnnotationDAG) {
        /**
         * Computes the information content of a set of nodes, returning it in a
         * map. Computed entries are cached for future usages.
         */
        Map<GOTerm, Float> informationContent = new HashMap<GOTerm, Float>();
        float invNumAnnotationDAG = 1.0f / (float) numAnnotationDAG;
        for (GOTerm term : CA) {
            if (this.icCache.containsKey(term)) {
                informationContent.put(term, this.icCache.get(term));
            } else {
                float ic = annotations.countNumberOfGenesForGOTerm(term.getGOid()) * invNumAnnotationDAG;
                informationContent.put(term, ic);
                this.icCache.put(term, ic);
            }
        }
        return informationContent;
    }

    private List<GOTerm> sortByDescIC(Map<GOTerm, Float> values) {
        List<TermWithIC> l = new ArrayList<TermWithIC>();
        for (Map.Entry<GOTerm, Float> e : values.entrySet()) {
            l.add(new TermWithIC(e.getKey(), e.getValue()));
        }
        Collections.sort(l, Collections.reverseOrder());
        List<GOTerm> sortedTerms = new ArrayList<GOTerm>();
        for (TermWithIC t : l) {
            sortedTerms.add(t.term);
        }
        return sortedTerms;
    }

    //Finds the disjunctive common ancestors for the two terms 'target1' & 'target2' in the ontology specified by 'dag'
    private float shareGraSM(GOTerm c1, GOTerm c2, int dag) {
        //Common ancestors
        Set<GOTerm> Anc = commonAncestors(c1, c2);
        Set<GOTerm> CommonDisjAnc = new HashSet<GOTerm>();
        Map<GOTerm, Float> IC = computeIC(Anc, this.maxAnnotationNumber[dag]);
        List<GOTerm> ancestorsSortedByICDesc = this.sortByDescIC(IC);

        for (GOTerm a : ancestorsSortedByICDesc) {
            boolean isDisj = true;

            for (GOTerm cda : CommonDisjAnc) {
                isDisj = isDisj && (DisjAnc(c1, a, cda) || DisjAnc(c2, a, cda));
                if (!isDisj) {
                    break;
                }
            }

            if (isDisj) {
                CommonDisjAnc.add(a);
            }
        }

        float shared = 0.0f;
        for (GOTerm cda : IC.keySet()) {
            shared += IC.get(cda);
        }
        return shared / (float) CommonDisjAnc.size();
    }

    @Override
    public Matrix calculateGeneWiseSemanticSimilarity(int ontology) throws IOException, OutOfMemoryError {
        return super.geneWiseSimilarityByMaximum(ontology);
    }

    @Override
    public Matrix calculateTermWiseSemanticSimilarity(int ontology) throws IOException, OutOfMemoryError {
        assert (ontology >= 0 && ontology < 3);
        final int N = numGOtermsPerOntology[ontology];
        Matrix result = new Matrix(N, N);

        for (int i = 0; i < N - 1; i++) {
            for (int j = i; j < N; j++) { //Semantic similarity calculated based upon disjunctive common ancestors
                float val = shareGraSM(this.matrixAxis[ontology][i], this.matrixAxis[ontology][j], ontology);
                result.set(i, j, val);
                result.set(j, i, val);
            }
        }
        logwriter.showMessage("Completed HSM for " + shortOntologyName[ontology]);
        return result;
    }

}
