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
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import util.TinyLogger;

/**
 *
 * @author Samuel Heron
 */
//Implements the simGraSM semantic similarity measure
public class simGraSM extends HSM {

    public simGraSM(GOTerm[] allTerms, String[] genes, String[][] goIds, GOTerm[][] axis, String[] targets, Assignment annotations, String[] relations, TinyLogger logw) {
        super(allTerms, genes, goIds, axis, targets, annotations, relations, logw);
    }

    //Returns the common ancestors of the GO terms 'target1' & 'target2'
    private HashSet<GOTerm> commonAncestors(GOTerm target1, GOTerm target2) {
        HashSet<GOTerm> CA = new HashSet<GOTerm>(target1.getAncestors()); //Common Ancestors
        CA.retainAll(target2.getAncestors());
        return CA;
    }

    //Finds the disjunctive common ancestors for the two terms 'target1' & 'target2' in the ontology specified by 'dag'
    private double disjunctiveCommonAncestors(GOTerm target1, GOTerm target2, int dag) {
        HashSet<GOTerm> DCA = new HashSet<GOTerm>(); //Disjunctive common ancestors
        HashSet<GOTerm> CA = commonAncestors(target1, target2); //Common ancestors
        HashSet<GOTerm> DAunion = new HashSet<GOTerm>(disjunctiveAncestors(target1)); //Union of disjunctive ancestors
        DAunion.addAll(disjunctiveAncestors(target2));

        //get DCA's
        // Common ancestors - DAunion
        CA.removeAll(DAunion);
        List<GOTerm> CAA = new ArrayList<GOTerm>(CA);
        double ic[] = new double[CAA.size()];
        for (int i = 0; i < CAA.size(); ++i) {
            ic[i] = annotations.countNumberOfGenesForGOTerm(CAA.get(i).getGOid()) / this.maxAnnotationNumber[dag];
        }

        for (int i = 0; i < CAA.size() - 1; ++i) {
            for (int j = i + 1; j < CAA.size(); ++j) {
                if (ic[i] <= ic[j]) {
                    DCA.add(CAA.get(i));
                    DCA.add(CAA.get(j));
                }
            }
        }

//get similarity value
        if (DCA.isEmpty()) //e.g root & root, no DA's, only CA = root etc.
        {
            return 0;
        } else {
            int counter = 0; //counts number of disjunctive common ancestors
            double semSimVal = 0.0;
            for (GOTerm dca : DCA) {
                semSimVal += -Math.log(annotations.countNumberOfGenesForGOTerm(dca.getGOid()) / this.maxAnnotationNumber[dag]); //Entropy value
                counter++;
            }
            semSimVal = semSimVal / counter;
            return semSimVal;
        }
    }
//Returns the disjunctive ancestors of a GO term; 'target'

    private HashSet<GOTerm> disjunctiveAncestors(GOTerm target) {
        HashSet<GOTerm> DA = new HashSet<GOTerm>(); //Disjunctive ancestors
        Set<GOTerm> ancestors = target.getAncestors();
        GOTerm[] ancs = ancestors.toArray(new GOTerm[0]);
        GOTerm a1, a2;
        for (int i = 0; i < ancs.length; i++) //For each ancestor
        {
            a1 = ancs[i];
            boolean conditionsMet1 = false;
            boolean conditionsMet2 = false;
            for (int j = i; j < ancs.length; j++) //Compare it against each other ancestor. No need to repeat results so j starts at i
            {
                a2 = ancs[j];
                for (HashSet<GOTerm> path : getPaths(a1, target)) //Fetches the paths between the ancestor and the target GO term
                {
                    if (path.contains(a2) == false) //If there is a path that doesn't contain the other ancestor then the first condition is met
                    {
                        conditionsMet1 = true;
                        break;
                    }
                }
                for (HashSet<GOTerm> path : getPaths(a2, target)) //Fetches the paths between the ancestor and the target GO term
                {
                    if (path.contains(a1) == false) //If there is a path that doesn't contain the other ancestor then the second condition is met
                    {
                        conditionsMet2 = true;
                        break;
                    }
                }
                if (conditionsMet1 == true && conditionsMet2 == true) //If both conditions are met then they are disjunctive ancestors
                {
                    DA.add(a1);
                    DA.add(a2);
                }
            }
        }
        return DA;
    }

    //Recursive method to find paths from a given GO term 'target' to another 'ancestor'. 'path' is the current path, 'allPaths' are all found paths between the terms
    private HashSet<HashSet<GOTerm>> pathFinder(GOTerm ancestor, GOTerm target, HashSet<GOTerm> path, HashSet<HashSet<GOTerm>> allPaths) {
        for (String rel : this.relations) {
            for (GOTerm parent : target.getParentsForRelation(rel)) {
                if (target.getParentsForRelation(rel).contains(ancestor) == true) //If a path is found from target to ancestor, add it to the path list
                {
                    allPaths.add(path);
                    return allPaths;
                } else if (target.isRoot(rel) == true) //If the current node is the root of an ontology, a path has not been found
                {
                    return allPaths;
                } else //Otherwise carry on looking
                {
                    path.add(parent);
                    pathFinder(ancestor, parent, path, allPaths);
                }
            }
        }
        return allPaths; //Return the paths found
    }

    //Returns the paths within the ontology between a given term 'target' and its ancestor; 'ancestor'
    private HashSet<HashSet<GOTerm>> getPaths(GOTerm ancestor, GOTerm target) {
        HashSet<HashSet<GOTerm>> resultantPaths = new HashSet<HashSet<GOTerm>>();
        HashSet<GOTerm> path = new HashSet<GOTerm>();
        resultantPaths = pathFinder(ancestor, target, path, resultantPaths);
        return resultantPaths;
    }

    @Override
    public RealMatrix calculateGeneWiseSemanticSimilarity(int ontology) throws IOException, OutOfMemoryError {
        return super.geneWiseSimilarityByMaximum(ontology);
    }

    @Override
    public RealMatrix calculateTermWiseSemanticSimilarity(int ontology) throws IOException, OutOfMemoryError {
        assert (ontology >= 0 && ontology < 3);
        final int N = numGOtermsPerOntology[ontology];
        RealMatrix result = new Array2DRowRealMatrix(N, N);

        for (int i = 0; i < N - 1; i++) {
            for (int j = i; j < N; j++) { //Semantic similarity calculated based upon disjunctive common ancestors
                double val = disjunctiveCommonAncestors(this.matrixAxis[ontology][i], this.matrixAxis[ontology][j], ontology);
                result.setEntry(i, j, val);
                result.setEntry(j, i, val);
            }
        }
        logwriter.showMessage("Completed HSM for " + shortOntologyName[ontology]);
        return result;
    }
}
