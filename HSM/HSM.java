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
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import util.TinyLogger;

/**
 *
 * @author Samuel Heron
 */
//This class is merely an interface for the HSM methods through which common functionality is shared
public abstract class HSM {

    //protected OpenMapRealMatrix[] GOadjacencyMatrices; //0=bp, 1=mf, 2=cc
    protected GOTerm[][] matrixAxis; //0=bp, 1=mf, 2=cc
    protected Assignment annotations; //annotation and obo details
    protected double[] maxAnnotationNumber; //the maximum number of annotations on any node
    //variables for genewise:
    protected String[] genes; //All gene IDs
    protected String[][] goIdsByGene; //A listing of GO terms by the Genes that are mapped to them
    protected GOTerm[] allTerms; //stores all GOTerms from all ontologies so that it can be searched to identify GOTerms
    protected HashMap<String, Integer> ontologyFromGOTerm; //a set of GO IDs truncated to their ontology id for easy identification when working with genes (as the GO terms will not be of GOTerm type)
    protected HashMap<String, Integer> indexFromGOTerm; //a set of GO IDs truncated to their ontology id for easy identification when working with genes (as the GO terms will not be of GOTerm type)
    protected HashMap<String, GOTerm> goTermFromID;
    protected String[] relations;
    protected TinyLogger logwriter;//Used for writing messages to the log file
    protected int numGOtermsPerOntology[];
    protected final int BIOLOGICAL_PROCESS = 0;
    protected final int MOLECULAR_FUNCTION = 1;
    protected final int CELLULAR_COMPONENT = 2;
    protected final static String[] shortOntologyName = {"BP", "MF", "CC"};
    protected final static String[] longOntologyName = {"Biological Process", "Molecular Function", "Cellular Component"};
    private String[] computedGenes;
    protected boolean isAGraphBasedMeasure;

    // Constructor for HSM genewise, takes a listing of all GO terms, a listing of all genes, a mapping 
    // of gene IDs to GO terms, whether or not a log file is to be written, the axis for the adjacencies, any specific terms, 
    // the adjacency matrices, the annotations, a choice of ontology & the logfile writer as parameters 
    // ArrayList<GOTerm> interestingTerms
    protected HSM(GOTerm[] allTerms, String[] genes, String[][] goIds, GOTerm[][] axis, String[] targets, /*OpenMapRealMatrix[] GO,*/ Assignment annotations, String[] relations, TinyLogger logw) {
        //initialisation of variables

        this.genes = genes;
        this.goIdsByGene = goIds;
        this.allTerms = allTerms;
        this.matrixAxis = axis;
        this.annotations = annotations;
        this.logwriter = logw;
        this.relations = relations;

        if (targets != null) {
            stripDownGeneData(targets);
        }

        numGOtermsPerOntology = new int[3];

        // ontology from GO Term
        this.ontologyFromGOTerm = new HashMap<String, Integer>();
        this.goTermFromID = new HashMap<String, GOTerm>();
        for (GOTerm go : this.allTerms) {
            this.goTermFromID.put(go.getGOid(), go);
            if (go.getOntology().getName().equals("biological_process")) {
                ontologyFromGOTerm.put(go.getGOid(), BIOLOGICAL_PROCESS);
                numGOtermsPerOntology[BIOLOGICAL_PROCESS]++;
            } else if (go.getOntology().getName().equals("molecular_function")) {
                ontologyFromGOTerm.put(go.getGOid(), MOLECULAR_FUNCTION);
                numGOtermsPerOntology[MOLECULAR_FUNCTION]++;
            } else {
                ontologyFromGOTerm.put(go.getGOid(), CELLULAR_COMPONENT);
                numGOtermsPerOntology[CELLULAR_COMPONENT]++;
            }
        }

        this.maxAnnotationNumber = getMaxAnnotations();

        // index from GO term
        this.indexFromGOTerm = new HashMap<String, Integer>();
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < matrixAxis[i].length; j++) {
                GOTerm term = matrixAxis[i][j];
                this.indexFromGOTerm.put(term.getGOid(), j);
            }
        }

        isAGraphBasedMeasure = false; // not a graph-based measure by default
    }

    //Returns the largest annotation value for normalisation purposes
    private double[] getMaxAnnotations() {

        double[] maxAnnoNo = new double[]{0, 0, 0};
        for (int i = 0; i < 3; ++i) {
            for (String term : annotations.getColumnIdentifiers()) {
                if (this.ontologyFromGOTerm.containsKey(term) && this.ontologyFromGOTerm.get(term) == i) {
                    maxAnnoNo[i] = Math.max(maxAnnoNo[i], annotations.countNumberOfGenesForGOTerm(term));
                }
            }
        }
        return maxAnnoNo;
    }

    public abstract Matrix calculateGeneWiseSemanticSimilarity(int ontology) throws IOException, OutOfMemoryError;

    public abstract Matrix calculateTermWiseSemanticSimilarity(int ontology) throws IOException, OutOfMemoryError;

    //Determine the lowest common ancestor of two terms & normalise the probability by the largest annotation value in the relevant ontology
    protected double lowestCommonAncestor(Set<GOTerm> ancestorsOne, Set<GOTerm> ancestorsTwo, int dag) {
        List<GOTerm> commonAncestors = new ArrayList<GOTerm>();
        int annoCount = Integer.MAX_VALUE;
        //get common ancestors
        for (GOTerm ancestor : ancestorsOne) {
            if (ancestorsTwo.contains(ancestor)) {
                commonAncestors.add(ancestor);
            }
        }
        GOTerm LCA = null;
        //find ancestor with lowest information content by referencing the common ancestors to the annotation assignment
        for (GOTerm comAnc : commonAncestors) {
            //if ancestor has less annotations than the previous lowest, store it as the LCA
            final int cnt = annotations.countNumberOfGenesForGOTerm(comAnc.getGOid());
            if (cnt < annoCount || LCA == null) {
                annoCount = cnt;
                LCA = comAnc;
            }
        }
        return (double) annoCount / maxAnnotationNumber[dag];
    }

    public int getNumGOTermsPerOntology(int ontology) {
        return this.numGOtermsPerOntology[ontology];
    }

    //returns a GOTerm's index from a given String GO ID
    protected int getOntologyFromGOTerm(String id) {
        if (this.ontologyFromGOTerm.containsKey(id)) {
            return this.ontologyFromGOTerm.get(id);
        } else {
            return -1;
        }
    }

    protected int getGOTermIndex(String id) {
        return this.indexFromGOTerm.get(id);
    }

    public String[] getSubSetGenes() {
        return this.computedGenes;
    }

    // removes unused gene2GO data & replace gene listing with the target genes
    private void stripDownGeneData(String[] targetGenes) {
        String[][] newGene2GO = new String[targetGenes.length][];
        int index = 0;
        for (String target : targetGenes) {
            for (int i = 0; i < this.genes.length; i++) {
                if (target.equals(this.genes[i]) == true) {
                    newGene2GO[index] = this.goIdsByGene[i];
                    index++;
                }
            }
        }
        this.genes = targetGenes;
        this.goIdsByGene = newGene2GO;
        boolean nulltest = false;
        for (int i = 0; i < newGene2GO.length; i++) {
            if (newGene2GO[i] == null) {
                nulltest = true;
                break;
            }
        }
        if (nulltest == true) {
            try {
                this.logwriter.logAndCloseWriter("############ERROR: Gene ids not found in GOA file");
                System.err.println("ERROR: Gene IDs entered could not be found in the GOA file supplied");
                System.exit(-1);
            } catch (IOException io) {
                System.err.println("ERROR: Gene IDs  entered could not be found in the GOA file supplied");
                System.exit(-1);
            }
        }
    }

    /**
     * This method tells if the HSM is a graph-based similarity measure
     * (Pesquita et al, 2008). By default this is false (the variable
     * "isAGraphBasedMeasure" should change its value to true for those cases in
     * which we are dealing with such a measure (e.g. simGIC).
     *
     * @returns true if the measure is a graph-based one, false otherwise
     */
    public boolean isAGraphBasedMeasure() {
        return isAGraphBasedMeasure;
    }

    protected float matrixMax(final Matrix mat) {
        float maxi = Float.NEGATIVE_INFINITY;
        final int rows = mat.getRowDimension();
        final int cols = mat.getColumnDimension();

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                float v = mat.get(i, j);
                maxi = maxi < v ? v : maxi;
            }
        }
        return maxi;
    }

    protected Matrix geneWiseSimilarityByMaximum(int ontology) throws IOException, OutOfMemoryError {
        //compute the semantic similarity
        Matrix termWise = this.calculateTermWiseSemanticSimilarity(ontology);

        Map<Integer, Set<Integer>> dominates = this.computeDominancies(termWise);

        Map<String, int[]> goIdsPerGene = new HashMap<String, int[]>();

        this.logwriter.showMessage("Extracting GO term ids associated to each gene... ");
        int all = 0, removed = 0;
        
        for (String gene : this.annotations.getRowIdentifiers()) {

            Set<Integer> ids = new HashSet<Integer>();
            Set<Integer> blacklisted = new HashSet<Integer>();

            for (String goTerm : this.annotations.getGOTermScoresForProteinId(gene).keySet()) {

                if (this.ontologyFromGOTerm.containsKey(goTerm) && this.ontologyFromGOTerm.get(goTerm) == ontology) {
                    //goterm index in HSM matrix.
                    int id = this.indexFromGOTerm.get(goTerm);

                    for (int id_in : ids) {
                        if (dominates.get(id).contains(id_in)) {
                            blacklisted.add(id_in);
                        } else if (dominates.get(id_in).contains(id)) {
                            blacklisted.add(id);
                            break;
                        }
                    }

                    ids.add(id);
                }
            }
            
            all+= ids.size();
            removed += blacklisted.size();

            if (!ids.isEmpty()) {

                int arrayIds[] = new int[ids.size() - blacklisted.size()];
                int i = 0;
                for (int val : ids) {
                    if (!blacklisted.contains(val)) {
                        arrayIds[i] = val;
                        i++;
                    }
                }

                Arrays.sort(arrayIds);
                goIdsPerGene.put(gene, arrayIds);
            }
        }
        
        double perc = (double)removed / (double) all * 100.0;
        
        this.logwriter.showMessage("Removed " + perc + "% of all associations because of dominancy");

        //filter out all the genes annotated to the desired ontology
        ArrayList<String> selectedGenes = new ArrayList<String>(goIdsPerGene.keySet());
        //sort the genes in alphabetical order.
        Collections.sort(selectedGenes);

        final int NUM_GENES_ONTOLOGY = selectedGenes.size();

        this.computedGenes = new String[NUM_GENES_ONTOLOGY];
        selectedGenes.toArray(this.computedGenes);
        
        this.logwriter.showMessage("Computing genewise semantic similarity by maximum (" + NUM_GENES_ONTOLOGY + " genes)");

        Matrix result = new Matrix(NUM_GENES_ONTOLOGY, NUM_GENES_ONTOLOGY);

        //which pair of terms annoating the genes is the most similar

        for (int i = 0; i < NUM_GENES_ONTOLOGY; ++i) {
            //get genes annotating the first gene
            final int[] goTerms_i = goIdsPerGene.get(selectedGenes.get(i));

            for (int j = i; j < NUM_GENES_ONTOLOGY; ++j) {
                //get genes annotating the second gene
                final int[] goTerms_j = goIdsPerGene.get(selectedGenes.get(j));

                float max = this.matrixMax(termWise.getMatrix(goTerms_i, goTerms_j));

                // set matrix values
                result.set(i, j, max);
                result.set(j, i, max);
            }
        }

        logwriter.log("Completed HSM for " + shortOntologyName[ontology]);
        System.out.println("Completed HSM for Ontology: " + longOntologyName[ontology]);

        return result;
    }

    private Map<Integer, Set<Integer>> computeDominancies(Matrix x) {
        final int m = x.getRowDimension();
        final int n = x.getRowDimension();

        Map<Integer, Set<Integer>> dominates = new HashMap<Integer, Set<Integer>>();

        for (int i = 0; i < m; ++i) {
            if (!dominates.containsKey(i)) {
                dominates.put(i, new HashSet<Integer>());
            }

            for (int j = i + 1; j < m; ++j) {
                if (!dominates.containsKey(j)) {
                    dominates.put(j, new HashSet<Integer>());
                }

                boolean greater = true, lesser = true;

                for (int k = 0; k < n; ++k) {
                    greater = greater && x.get(i, k) >= x.get(j, k);
                    lesser = lesser && x.get(i, k) < x.get(j, k);

                    if (!greater && !lesser) {
                        break;
                    }
                }

                if (greater) {
                    dominates.get(i).add(j);
                }
                if (lesser) {
                    dominates.get(j).add(i);
                }
            }
        }

        /*
        int k = 0;
        for (int i = 0; i < m; ++i) {
            if (dominates.get(i).size() > 10) {
                ++k;
                System.err.println("Row " + i + " dominates " + dominates.get(i).size());
            }
        } */
        return dominates;
    }
}
