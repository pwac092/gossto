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
package ISM;

import GOtree.AnnotationFile;
import GOtree.Assignment;
import GOtree.GOTerm;
import GOtree.GeneOntology;
import GOtree.GeneOntologyParser;
import GOtree.Ontology;
import GOtree.Propagation;
import GOtree.PropagationStrategies.PropagationMaxWithinTrees;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.apache.commons.math3.linear.OpenMapRealMatrix;
import util.TinyLogger;

/**
 *
 * @author Samuel Heron
 */
//This class deals with preparing and importing all relevant source data, it provides a buffer between the main class and the classes dealing
//solely with the source data
public class GOtreeInterfacer {

    private GeneOntology GO; //the gene ontology itself
    private OpenMapRealMatrix BPAdjacencyMatrix, MFAdjacencyMatrix, CCAdjacencyMatrix; //the adjacency matrices for the three ontologies
    private String[] GO_relations, geneIDs; //The relations linking the GO terms that we will be using
    private Assignment annotations; //the mapping of the annotations to the GO terms
    private GOTerm[] bpAxisOrd, mfAxisOrd, ccAxisOrd; //The axis for the adjacency matrices
    private String[][] goIDsByGene; //Listing of GO IDs that a gene is annotated to
    private HashMap<GOTerm, Double> weightedAnnotations; //Contains the weighted annotation mappings used for calculating HSM results to be used for calculating ISM results
    private TinyLogger logwriter;//Used for writing messages to the log file
    Set<String> ontologiesToProcess;
    private HashMap<Integer, Integer> numAnnotations = new HashMap<Integer, Integer>();

    //All methods are called from this constructor as no user interaction is required.
    //'OBOpath' contains the file path for the OBO file, 'annoPath' contains the file path for the GOA file, 'relations' contains the GO relations we are using, 'evidenceCodes' contains
    //the evidence codes we're using, 'strategyChoice' is the propagation strategy to use, however this is always 1, and the two variables are
    //used for writing the log file if necessary & are used for thus as parameters in the other methods within this class
    public GOtreeInterfacer(String OBOpath, String annoPath, String[] relations, String[] evidenceCodes, int strategyChoice, String dagChoice, TinyLogger logw) throws FileNotFoundException, IOException {
        //all methods called from the constructor, no need to reference anything but the results produced and made available by the 'getter' methods
        //instantiates the log variables if necessary
        logwriter = logw;
        this.ontologiesToProcess = new HashSet<String>();

        if (dagChoice.equalsIgnoreCase("ALL")) {
            ontologiesToProcess.add("BP");
            ontologiesToProcess.add("MF");
            ontologiesToProcess.add("CC");
        } else if (dagChoice.equalsIgnoreCase("BP")) {
            ontologiesToProcess.add("BP");
        } else if (dagChoice.equalsIgnoreCase("MF")) {
            ontologiesToProcess.add("MF");
        } else if (dagChoice.equalsIgnoreCase("CC")) {
            ontologiesToProcess.add("CC");
        }
        this.GO_relations = relations;
        this.weightedAnnotations = new HashMap<GOTerm, Double>();
        logwriter.log("Fetching GO");
        generate_GO(OBOpath); //Parse the OBO file
        logwriter.log("GO parsed");
        logwriter.log("Loading annotation");
        load_Annotation(annoPath, evidenceCodes); //Parse the GOA file
        logwriter.log("Annotation loaded");
        logwriter.log("Fetching Gene to GO mappings");
        getGOtermsForProteins(); //Get gene to GO term mappings
        logwriter.log("Creating adjacency matrices");
        try {
            generate_AdjacencyMatrices(this.GO); //Create adjacency matrices
            logwriter.log("Propagating annotations");
            propagator(strategyChoice); //Up propagate annotations
            logwriter.log("Propagation complete");
            logwriter.log("Stripping data");
            stripper(); //strip down adjacency matrices to valuable data & likewise for their matrices
        } catch (java.lang.OutOfMemoryError oome) {
            this.logwriter.logAndCloseWriter("############## ERROR: Out of memory Error of type: " + oome.getMessage());
            System.err.println("ERROR: Java has run out of memory. Memory Type: " + oome.getMessage());
            System.exit(-1);
        }
    }
    
    public boolean annotationFileIsEmpty()
    {
        return this.annotations.sizeTerms() == 0;
    }

    //parses the OBO file specified by 'OBOpath'
    private void generate_GO(String OBOpath) throws FileNotFoundException, IOException {
        GeneOntologyParser parser = new GeneOntologyParser();
        GOTerm.setRelations(this.GO_relations);
        GeneOntology ontology = parser.readFromOBOFile(OBOpath, false);
        this.GO = ontology;
    }

    //parses the GOA file specified by 'annoPath' & uses the evidence Codes specified by 'evidenceCodes'
    private void load_Annotation(String annoPath, String[] evidenceCodes) throws FileNotFoundException, IOException {
        AnnotationFile annoFile = new AnnotationFile();
        this.annotations = annoFile.readAnnotationFile(annoPath, evidenceCodes);
    }

    //Up propagates the annotations according to the propagation strategy specified by 'strategyChoice' (always PropagationMaxWithinTrees)
    private void propagate_Annotations(int strategyChoice) {
        Propagation propagator = new Propagation();
        PropagationMaxWithinTrees strategy = new PropagationMaxWithinTrees();
        Assignment propagated_annotations = propagator.propagateAssignment(annotations, GO, strategy);
        this.annotations = propagated_annotations;
    }

    class CustomComparator implements Comparator<GOTerm> {

        @Override
        public int compare(GOTerm t, GOTerm t1) {
            return t.getNumericId() - t1.getNumericId();
        }
    }

    //order the GO terms, specified by 'terms', in an ascending manner by their numerical id
    private GOTerm[] orderGOTerms(List<GOTerm> terms) {
        GOTerm[] termArray = new GOTerm[terms.size()];
        terms.toArray(termArray);
        Arrays.sort(termArray);
        return termArray;
    }
//Creates the adjacency matrices that represent the structure of the Gene Ontologies specified by 'go'

    private void generate_AdjacencyMatrices(GeneOntology go) throws IOException {
        //fetch ontologies
        Ontology bp, mf, cc;
        this.bpAxisOrd = new GOTerm[0];
        this.mfAxisOrd = new GOTerm[0];
        this.ccAxisOrd = new GOTerm[0];

        Map<Integer, Integer> indexBP = new HashMap<Integer, Integer>();
        Map<Integer, Integer> indexMF = new HashMap<Integer, Integer>();
        Map<Integer, Integer> indexCC = new HashMap<Integer, Integer>();


        for (String ontology : this.ontologiesToProcess) {
            if (ontology.trim().equalsIgnoreCase("BP")) {

                bp = go.getOntology("biological_process");
                List<GOTerm> bpAxis = bp.getSetOfTerms();
                this.bpAxisOrd = new GOTerm[bpAxis.size()];
                this.bpAxisOrd = orderGOTerms(bpAxis);

                int i = 0;
                for (GOTerm term : bpAxisOrd) {
                    indexBP.put(term.getNumericId(), i);
                    i++;
                }

                this.BPAdjacencyMatrix = new OpenMapRealMatrix(bp.size(), bp.size());

                //set adjacencies' values ontology by ontology
                logwriter.showMessage("Creating BP adjacency");

                int currentIndex = 0;
                for (GOTerm term : this.bpAxisOrd) { //for all the goterms
                    for (String rel : GO_relations) { //check for all relations
                        for (GOTerm parent : term.getParentsForRelation(rel)) { //get all the parents
                            //int idx = findIndex(this.bpAxisOrd, parent); //find the corresponding index for the goterm
                            //this checks for parents that might belong to a different ontology

                            //if (idx != -1) {

                            if (indexBP.containsKey(parent.getNumericId())) {
                                int idx = indexBP.get(parent.getNumericId());
                                this.BPAdjacencyMatrix.setEntry(currentIndex, idx, 1); //set a one indicating the connection in the 
                            }

                            //}
                        }
                    }
                    currentIndex++;
                }

            } else if (ontology.trim().equalsIgnoreCase("MF")) {
                mf = go.getOntology("molecular_function");
                List<GOTerm> mfAxis = mf.getSetOfTerms();

                this.mfAxisOrd = new GOTerm[mfAxis.size()];
                this.mfAxisOrd = orderGOTerms(mfAxis);

                int i = 0;
                for (GOTerm term : mfAxisOrd) {

                    indexMF.put(term.getNumericId(), i);
                    i++;
                }

                this.MFAdjacencyMatrix = new OpenMapRealMatrix(mf.size(), mf.size());

                //set adjacencies' values ontology by ontology
                logwriter.showMessage("Creating MF adjacency");

                int currentIndex = 0;
                for (GOTerm term : this.mfAxisOrd) {
                    for (String rel : GO_relations) {
                        for (GOTerm parent : term.getParentsForRelation(rel)) {
                            // int idx = findIndex(this.mfAxisOrd, parent);
                            //this checks for parents that might belong to a different ontology
                            //if (idx != -1) {
                            //    this.MFAdjacencyMatrix.setEntry(currentIndex, idx, 1);
                            //}
                            if (indexMF.containsKey(parent.getNumericId())) {
                                int idx = indexMF.get(parent.getNumericId());
                                this.MFAdjacencyMatrix.setEntry(currentIndex, idx, 1); //set a one indicating the connection in the 
                            }
                        }
                    }
                    currentIndex++;
                }

            } else if (ontology.trim().equalsIgnoreCase("CC")) {
                cc = go.getOntology("cellular_component");
                List<GOTerm> ccAxis = cc.getSetOfTerms();
                this.ccAxisOrd = new GOTerm[ccAxis.size()];
                this.ccAxisOrd = orderGOTerms(ccAxis);

                int i = 0;
                for (GOTerm term : ccAxisOrd) {
                    indexCC.put(term.getNumericId(), i);
                    i++;
                }

                this.CCAdjacencyMatrix = new OpenMapRealMatrix(cc.size(), cc.size());
                //set adjacencies' values ontology by ontology
                logwriter.showMessage("Creating CC adjacency");

                int currentIndex = 0;
                for (GOTerm term : this.ccAxisOrd) {
                    for (String rel : GO_relations) {
                        for (GOTerm parent : term.getParentsForRelation(rel)) {
                            //int idx = findIndex(this.ccAxisOrd, parent);
                            //this checks for parents that might belong to a different ontology
                            //if (idx != -1) {
                            //    this.CCAdjacencyMatrix.setEntry(currentIndex, idx, 1);
                            //}
                            if (indexCC.containsKey(parent.getNumericId())) {
                                int idx = indexCC.get(parent.getNumericId());
                                this.CCAdjacencyMatrix.setEntry(currentIndex, idx, 1); //set a one indicating the connection in the 
                            }
                        }
                    }
                    currentIndex++;
                }

            }
        }
    }

    //Either calls the propagate_Annotations() method if the end product is an HSM or uses the weighted propagation mechanism within this method
    //if the end product is an ISM
    //'ISM' determines the propagation method & 'strategyChoice' determines the method of HSM propagation (always going to be propagationMax)
    private void propagator(int strategyChoice) {
        System.out.println("Propagating");
        propagate_Annotations(strategyChoice);
    }
//Strips down adjacencies and their axis to only relevant information, the end product of this execution (HSM or ISM) specified by 'ISM' determines
//how the stripping will be done

    private void stripper() throws IOException {
        logwriter.log("Stripping adjacencies");
        //stripping for HSM
        //strip down adjacencies and axis to those with annotations
        logwriter.log("Stripping axis");

        for (String tree : this.ontologiesToProcess) {
            if (tree.equalsIgnoreCase("BP")) {
                this.BPAdjacencyMatrix = stripDownAdjacency(this.BPAdjacencyMatrix, this.bpAxisOrd);
                this.bpAxisOrd = stripDownAxis(this.bpAxisOrd);
            } else if (tree.equalsIgnoreCase("MF")) {
                this.MFAdjacencyMatrix = stripDownAdjacency(this.MFAdjacencyMatrix, this.mfAxisOrd);
                this.mfAxisOrd = stripDownAxis(this.mfAxisOrd);
            } else if (tree.equalsIgnoreCase("CC")) {
                this.CCAdjacencyMatrix = stripDownAdjacency(this.CCAdjacencyMatrix, this.ccAxisOrd);
                this.ccAxisOrd = stripDownAxis(this.ccAxisOrd);
            }
        }

    }

    private int getNumAnnotationsForGOTerm(GOTerm term) {
        final int numId = term.getNumericId();
        if (this.numAnnotations.containsKey(numId)) {
            return this.numAnnotations.get(numId);
        } else {
            int num = this.annotations.countNumberOfGenesForGOTerm(term.getGOid());
            this.numAnnotations.put(numId, num);
            return num;
        }
    }

    //remove columns and rows of terms without annotations for the matrix 'adjacency' with the axis 'axis'
    private OpenMapRealMatrix stripDownAdjacency(OpenMapRealMatrix adjacency, GOTerm[] axis) {

        //get number of terms with annotations
        final int n = adjacency.getRowDimension();
        Set<Integer> interesting = new HashSet<Integer>();

        for (int i = 0; i < n; i++) {
            //if (annotations.countNumberOfGenesForGOTerm(axis[i].getGOid()) != 0) {
            if (this.getNumAnnotationsForGOTerm(axis[i]) != 0) {
                //counter++;
                interesting.add(i);
            }
        }
        final int counter = interesting.size();

        if (counter > 0) {

            OpenMapRealMatrix stripped = new OpenMapRealMatrix(counter, counter);
            //strip matrix down
            int rowIndex = 0;
            //final int m = adjacency.getColumnDimension();
            for (int i : interesting) {
                int columnIndex = 0;
                //if (annotations.countNumberOfGenesForGOTerm(axis[i].getGOid()) != 0) {
                //if (this.getNumAnnotationsForGOTerm(axis[i]) != 0) {
                //for (int j = 0; j < m; j++) {
                for (int j : interesting) {
                    //if (this.getNumAnnotationsForGOTerm(axis[j]) != 0) {
                    //if (annotations.countNumberOfGenesForGOTerm(axis[j].getGOid()) != 0) {                
                    stripped.setEntry(rowIndex, columnIndex, adjacency.getEntry(i, j));
                    columnIndex++;
                    //}
                }
                rowIndex++;
                //}
            }
            return stripped;
        } else {
            return null;
        }
    }
    //remove terms from the axis that don't have annotations for the axis; 'axis'

    private GOTerm[] stripDownAxis(GOTerm[] axis) {
        GOTerm[] stripped;
        int counter = 0;
        //get number of terms with annotations
        for (int i = 0; i < axis.length; i++) {
            if (annotations.countNumberOfGenesForGOTerm(axis[i].getGOid()) != 0) {
                counter++;
            }
        }
        stripped = new GOTerm[counter];
        counter = 0;
        //remove the terms without annotations
        for (int i = 0; i < axis.length; i++) {
            if (annotations.countNumberOfGenesForGOTerm(axis[i].getGOid()) != 0) {
                stripped[counter] = axis[i];
                counter++;
            }
        }
        return stripped;
    }

    //fetches the geneIDs and the gene-to-goterm mappings
    private void getGOtermsForProteins() {
        String[] genes = new String[annotations.sizeGenes()];
        String[][] goIDs = new String[annotations.sizeGenes()][];
        Map<String, Double> temp;
        for (int i = 0; i < annotations.sizeGenes(); i++) {
            genes[i] = annotations.getGeneFromId(i);
            temp = annotations.getGOTermScoresForProteinId(genes[i]);
            goIDs[i] = mapToString(temp);
        }
        this.geneIDs = genes;
        this.goIDsByGene = goIDs;
    }

    //converts a mapping; 'map', to a string
    private String[] mapToString(Map<String, Double> map) {
        String raw = map.toString();
        String[] interim = raw.split(",");
        String[] goIds = new String[interim.length];
        int counter = 0;
        for (String rawTerm : interim) {
            goIds[counter] = rawTerm.substring(1, 11);
            counter++;
        }
        return goIds;
    }

    public GeneOntology getGO() {
        return this.GO;
    }

    public Assignment getResults() {
        return this.annotations;
    }

    public GOTerm[] getBPaxis() {
        return this.bpAxisOrd;
    }

    public GOTerm[] getMFaxis() {
        return this.mfAxisOrd;
    }

    public GOTerm[] getCCaxis() {
        return this.ccAxisOrd;
    }

    public OpenMapRealMatrix getBPAdjacencyMatrix() {
        return this.BPAdjacencyMatrix;
    }

    public OpenMapRealMatrix getCCAdjacencyMatrix() {
        return this.CCAdjacencyMatrix;
    }

    public OpenMapRealMatrix getMFAdjacencyMatrix() {
        return this.MFAdjacencyMatrix;
    }

    public String[] getGO_relations() {
        return this.GO_relations;
    }

    public String[] getGeneIDs() {
        return this.geneIDs;
    }

    public String[][] getGoIdsByGene() {
        return this.goIDsByGene;
    }

    public HashMap<GOTerm, Double> getWeightedAnnoMappings() {
        return this.weightedAnnotations;
    }
}
