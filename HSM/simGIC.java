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
package HSM;

import GOtree.Assignment;
import GOtree.GOTerm;
import java.io.IOException;
import java.util.HashSet;
import org.apache.commons.math3.linear.OpenMapRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import util.TinyLogger;

/**
 *
 * @author Samuel Heron
 */
//Implements the simGIC semantic similarity measure
//get max no. of annotations on any one gene
public class simGIC extends HSM {

    public simGIC(GOTerm[] allTerms, String[] genes, String[][] goIds, GOTerm[][] axis, String[] targets,  Assignment annotations, String[] relations, TinyLogger logw) {
        super(allTerms, genes, goIds, axis, targets, annotations, relations, logw);
    }

    private double getMaxAnnotations() {
        double maxAnnoNo = 0;
        for (String term : annotations.getColumnIdentifiers()) {
            if (annotations.countNumberOfGenesForGOTerm(term) > maxAnnoNo) //if a node has more, then make that the most
            {
                maxAnnoNo = (double) annotations.countNumberOfGenesForGOTerm(term);
            }
        }
        //}
        return maxAnnoNo;
    }

    //Fetches the GOTerm instance of a term from its ID
    /*
     private GOTerm getGOTerm(String id) {
     String bp, mf;
     bp = id + "biological_process";
     mf = id + "molecular_function";
     if (idsAndDags.contains(bp) == true) {
     for (int i = 0; i < matrixAxis[0].length; i++) {
     if (matrixAxis[0][i].getGOid().equals(id) == true) {
     return matrixAxis[0][i];
     }
     }
     } else if (idsAndDags.contains(mf) == true) {
     for (int i = 0; i < matrixAxis[1].length; i++) {
     if (matrixAxis[1][i].getGOid().equals(id) == true) {
     return matrixAxis[1][i];
     }
     }
     } else {
     for (int i = 0; i < matrixAxis[2].length; i++) {
     if (matrixAxis[2][i].getGOid().equals(id) == true) {
     return matrixAxis[2][i];
     }
     }
     }
     try {
     this.logwriter.logAndCloseWriter("#######ERROR: GO Term doesn't exist");
     } catch (IOException e) {
     e.printStackTrace(System.err);
     }
     System.err.println("ERROR: GO term doesn't exist");
     System.exit(-1);
     return new GOTerm(null, null);
     }*/
    @Override
    public RealMatrix calculateTermWiseSemanticSimilarity(int ontology) throws IOException, OutOfMemoryError {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public RealMatrix calculateGeneWiseSemanticSimilarity(int ontology) throws IOException, OutOfMemoryError {
        System.out.println("#####SimGIC HSM#####");
        final int N = numGOtermsPerOntology[ontology];
        RealMatrix result = new OpenMapRealMatrix(N, N);
        final int NUM_GENES = this.genes.length;

        // GeneGoups for any gene has all the ancestors of the GO terms annotated
        // to that gene
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
        double maxAnno = getMaxAnnotations();

        for (int i = 0; i < NUM_GENES-1; i++) {
            if (!GeneGroups[i].isEmpty()) {
                for (int j = i; j < NUM_GENES; j++) {
                    if (!GeneGroups[j].isEmpty()) {
                        double intersectionVal = 0;

                        for (GOTerm go : GeneGroups[i])//get values that are unique to one or other set
                        {
                            if (!GeneGroups[j].contains(go)) {
                                double ICtemp = this.annotations.countNumberOfGenesForGOTerm(go.getGOid());
                                intersectionVal += -Math.log(ICtemp / maxAnno);
                            }
                        }
                        for (GOTerm go : GeneGroups[j]) {
                            if (!GeneGroups[i].contains(go)) {
                                double ICtemp = this.annotations.countNumberOfGenesForGOTerm(go.getGOid());
                                intersectionVal += -Math.log(ICtemp / maxAnno);
                            }
                        }
                        HashSet<GOTerm> union = new HashSet<GOTerm>();
                        double unionVal = 0;
                        union.addAll(GeneGroups[i]); //get union of the two sets
                        union.addAll(GeneGroups[j]);
                        for (GOTerm go : union) {
                            //get IC of each term in the intersection
                            double ICtemp = this.annotations.countNumberOfGenesForGOTerm(go.getGOid());
                            unionVal += -Math.log(ICtemp / maxAnno);
                        }
                        double semSim = intersectionVal / unionVal; //calculate similarity value
                        result.setEntry(i, j, semSim);
                        result.setEntry(j, i, semSim);
                    }
                }
            }
        }

        logwriter.log("Completed HSM for " + shortOntologyName[ontology]);
        System.out.println("Completed simGIC for Ontology : Biological Process" + longOntologyName[ontology]);

        return result;
    }
}
