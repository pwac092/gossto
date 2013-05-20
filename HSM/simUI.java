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
import java.util.HashSet;
import util.TinyLogger;

/**
 *
 * @author Samuel Heron, Alfonso E. Romero
 * 
 */
//Implements the simUI semantic similarity measure
public class simUI extends HSM {

    public simUI(GOTerm[] allTerms, String[] genes, String[][] goIds, GOTerm[][] axis, String[] targets, Assignment annotations, String[] relations, TinyLogger logw) {
        super(allTerms, genes, goIds, axis, targets, annotations, relations, logw);
        isAGraphBasedMeasure = true;
    }

    @Override
    public Matrix calculateGeneWiseSemanticSimilarity(int ontology) throws IOException, OutOfMemoryError {
        System.out.println("#####SimUI HSM#####");
        final int N = numGOtermsPerOntology[ontology];
        Matrix result = new Matrix(N, N);
        final int NUM_GENES = this.genes.length;

        HashSet<GOTerm> GeneGroups[] = new HashSet[NUM_GENES];
        for (int g = 0; g < NUM_GENES; g++) {
            GeneGroups[g] = new HashSet<GOTerm>();
        }
        //get all GOTerms and ancestors for each gene
        logwriter.showMessage("Getting gene2go info for " + shortOntologyName[ontology]);

        for (int g = 0; g < this.genes.length; g++) {
            for (String go : this.goIdsByGene[g]) {
                if (ontology == super.getOntologyFromGOTerm(go)) {
                    GeneGroups[g].addAll(super.goTermFromID.get(go).getAncestors());
                }
            }
        }

        logwriter.log("Getting gene2go Completed ");
        System.out.println("Calculating Semantic Similarity");
        //calculate semantic similarity between each gene 
        HashSet<GOTerm> union = new HashSet<GOTerm>();
        for (int i = 0; i < GeneGroups.length; i++) {
            if (GeneGroups[i].isEmpty() == false) {
                for (int j = 0; j < GeneGroups.length; j++) {
                    if (GeneGroups[j].isEmpty() == false) {
                        int intersectionVal = 0;
                        union.clear();
                        for (GOTerm go : GeneGroups[i]) //get values that are unique to one or other set
                        {
                            if (GeneGroups[j].contains(go) == false) {
                                intersectionVal++;
                            }
                        }
                        for (GOTerm go : GeneGroups[j]) {
                            if (GeneGroups[i].contains(go) == false) {
                                intersectionVal++;
                            }
                        }
                        union.addAll(GeneGroups[i]); //get union of the two sets
                        union.addAll(GeneGroups[j]);
                        float unionVal = union.size();
                        float semSim = (float) intersectionVal / unionVal; //calculate similarity value
                        result.set(i, j, semSim);
                    }
                }
            }
        }
        logwriter.showMessage("Completed HSM for " + shortOntologyName[ontology]);
        return result;
    }

    @Override
    public Matrix calculateTermWiseSemanticSimilarity(int ontology) throws IOException, OutOfMemoryError {
        throw new UnsupportedOperationException("Not supported yet.");
    }

}
