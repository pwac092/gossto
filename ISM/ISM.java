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
package ISM;

import GOtree.AnnotationFile;
import GOtree.GOTerm;
import GOtree.GeneOntologyException;
import Jama.Matrix;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;
import util.TinyLogger;

/**
 *
 * @author Samuel Heron, Alfonso E. Romero
 */
//Contains the instructions for running GOssTo, remove this class to use GOssTo as a library
public class ISM {

    /**
     * Interface object interacts with the user
     */
    private TerminalInterface ui;
    /**
     * Name of the OBO file which contains the Gene Ontology used
     */
    private String oboFile;
    /**
     * Name of the GOA file (annotations) used to assign functions to genes
     */
    private String goaFile;
    /**
     * Name of the Host Similarity measure used
     */
    private String hsmChoice;
    /**
     * Names of the ontologies used to compute the HSMs (and possibly ISMs)
     */
    private String dagChoice;
    /**
     * Name of the file where the HSM is written
     */
    private String hsmFileName;
    /**
     * Name of the file where the ISM is written
     */
    private String ismFileName;
    /**
     * Relations in the Gene Ontology which are used to up-propagate
     */
    private String[] chosenRelations;
    /**
     * Gene identifiers
     */
    private String[] geneIDs;
    /**
     * specific set of GO terms that are to be used in the execution (strings)
     */
    private String[] goIDs;
    /**
     * List of evidence codes considered in the GOA file
     */
    private String[] evidenceCodes;
    /**
     * Tells whether the ISM is going to be computed or not
     */
    private boolean isIsmToBeComputed;
    /**
     * Tells whether the measure is computed termwise
     */
    private boolean termWise;
    /**
     *
     */
    private boolean weightedJaccard;
    /**
     * Flag that tells us if the printing should be done in matrix style (true),
     * or in triplet style (false)
     */
    private int matrixStyle;
    /**
     * Tells us if we are using UniProtKB accesion numbers (true) or gene names
     * to identify proteins
     */
    private boolean useUniProtIds;
    /**
     * Logger used to output messages
     */
    TinyLogger logger;
    /**
     * Stores the execution parameters to write them at the top of the output
     * files
     */
    private ArrayList<String> notes;
    public static final int MATRIX_STYLE = 0;
    public static final int TRIPLET_STYLE = 1;
    public static final int BOTH_FILES = 2;

    /**
     * Empty constructor
     * @throws java.io.FileNotFoundException
     */
    public ISM() throws FileNotFoundException {
        logger = new TinyLogger();
        IoValidation.setLogger(logger);
        this.ui = new TerminalInterface();
        ui.welcomer(); //Prints welcome note & GNU GPL info
    }

    /**
     * Validates parameters given to the ISM application, either in the terminal
     * or via the prompt
     *
     * @param args arguments given in the command line to the program
     */
    private void validateParameters(String args[]) throws FileNotFoundException, IOException {
        ParameterValidator paramValidator;

        if (args.length > 0) //If any console parameters entered:
        {
            paramValidator = new ConsoleParameterValidator(args);
            paramValidator.validate(this.logger);
            this.setParametersConsole(paramValidator);

        } else //Request Parameters using prompt system
        {
            paramValidator = new PromptParameterValidator(ui);
            paramValidator.validate(this.logger);
            this.setParametersPrompt(paramValidator);
        }
        AnnotationFile.useUniProtIds(this.useUniProtIds);
        this.logger.log("All parameters validated, except GO terms");
    }

    private void setParametersConsole(ParameterValidator validator) {
        oboFile = validator.getOboFile();
        goaFile = validator.getGoaFile();
        hsmChoice = validator.getHsmChoice();
        dagChoice = validator.getDagChoice();
        hsmFileName = validator.getHsmFileName();
        ismFileName = validator.getIsmFileName();
        chosenRelations = validator.getChosenRelations();
        geneIDs = validator.getGeneIDs();
        goIDs = validator.getGoIDs();
        evidenceCodes = validator.getEvidenceCodes();
        notes = validator.getNotes();
        this.isIsmToBeComputed = validator.isIsmChoice();
        this.termWise = validator.isTermWise();
        this.weightedJaccard = validator.isWeightedJaccard();
        this.useUniProtIds = validator.isUseUniProtIds();
        this.matrixStyle = validator.getMatrixStyle();
    }

    private void setParametersPrompt(ParameterValidator validator) {
        this.setParametersConsole(validator);
    }

    /**
     * Imports the GO tree from file
     *
     * @return a GOtree_interfacer object with all the GO tree data
     */
    private GOtreeInterfacer importGOTree() throws FileNotFoundException, IOException, GeneOntologyException {

        logger.showMessage("#####Importing GO & Annotation Data#####");
        int propagationStrategy = 1; //Choice of propagation strategy, 1 as default, never changed.
        GOtreeInterfacer gti = new GOtreeInterfacer(this.oboFile, this.goaFile, this.chosenRelations, this.evidenceCodes, propagationStrategy, this.dagChoice, this.logger);
        logger.log("GOtree_Interfacer instantiated & executed");
        return gti;
    }

    /**
     * Returns the matrix axis (a structure where the GO terms are organised by
     * ontology) from the GO tree interfacer
     *
     * @param gti GO tree data imported by the previous step
     * @return a matrix of GO terms (adjacencies) where the first dimension is
     * one of each three ontologies (BP, MF, CC) and the second the list of
     * terms
     */
    private static GOTerm[][] getMatrixAxis(GOtreeInterfacer gti) {
        GOTerm[][] matrixAxis = new GOTerm[3][];
        matrixAxis[0] = gti.getBPaxis(); //fetch adjacency matrix axis
        matrixAxis[1] = gti.getMFaxis();
        matrixAxis[2] = gti.getCCaxis();

        return matrixAxis;
    }

    /**
     * Validates the GO terms
     *
     * @param matrixAxis matrix axis imported by the previous step
     * @return list of interesting GO term object considered by the user which
     * are going to be used to reduce the dimensionality of the matrices (might
     * be empty)
     */
    private ArrayList<GOTerm> validateGOTerms(GOTerm[][] matrixAxis) throws IOException {

        // converts chosen this.goIDs (if any) from Strings to a list of GO terms 
        // and validate them simultaneously
        ArrayList<GOTerm> goIDsAsGOTerm = new ArrayList<GOTerm>();
        ArrayList<String> termsNotFound = new ArrayList<String>();
        if (this.goIDs != null) {
            for (String id : this.goIDs) {
                boolean termFound = false;
                SEARCHING:
                for (int m = 0; m < 3; m++) {
                    for (GOTerm item : matrixAxis[m]) {
                        if (id.equals(item.getGOid())) {
                            goIDsAsGOTerm.add(item);
                            termFound = true;
                            break SEARCHING;
                        }
                    }
                }
                if (!termFound) {
                    termsNotFound.add(id);
                }
            }
            if (!termsNotFound.isEmpty()) {
                int num = termsNotFound.size();
                if (num > 1) {
                    logger.showMessage("Warning: " + num + " of the specified GO terms could not be found.\nDetails:");
                } else {
                    logger.showMessage("Warning: one of the specified GO terms could not be found. Details:");
                }
                for (String term : termsNotFound) {
                    logger.showMessage("\t+ " + term);
                }
            }


            /*
             if (this.goIDs.length != goIDsAsGOTerm.size()) {
             logger.logAndCloseWriter("############ ERROR: Invalid GOTerms found");
             System.err.println("ERROR: One or more goterm IDs entered do not exist");
             System.exit(-1);
             }*/
        }

        logger.log("Terms validated");
        return goIDsAsGOTerm;
    }

    /**
     * Contains the main running instructions for the program, takes console
     * parameters as an array of strings.
     * @param args arguments of the program
     */
    public static void main(String[] args) {
        try {

            //#####Proper Running of Program#####
            ISM ism = new ISM();

            // 1.- the parameters are validated
            ism.validateParameters(args);

            // 2.- the GO tree is imported and validated wrt annotation files
            GOtreeInterfacer gti = ism.importGOTree();
            GOTerm[][] matrixAxis = getMatrixAxis(gti);
            ArrayList<GOTerm> goTerms = ism.validateGOTerms(matrixAxis);

            // 3.- the semantic similarities are computed and written to disk
            ism.computeAndWriteSemanticSimilarities(gti, matrixAxis, goTerms);

            // 4.- Say goodbye!
            ism.farewell();

        } catch (OutOfMemoryError ex) {
            System.err.println("ERROR: Insufficient memory to run GOSSTO with the chosen parameter set.");
            System.err.println("Please, launch the Java Virtual Machine with at least 2 GB of memory,");
            System.err.println("by setting 'java -Xmx2G ... '. Check your systems documentation for");
            System.err.println("specific options.");
            System.exit(-1);
        } catch (GeneOntologyException ex) {
            System.err.println("ERROR: problem with the Gene Ontology file.");
            System.err.println(ex.getMy_message());
            System.exit(-1);
        } catch (IOException ex) {
            ex.printStackTrace(System.err);
            System.exit(-1);
        }
    }

    private void computeAndWriteSemanticSimilarities(GOtreeInterfacer gti, GOTerm[][] matrixAxis, ArrayList<GOTerm> goIDsAsGOTerm) throws IOException {

        // 1.- we compute the adjacency matrices for each Ontology.
        // the "adjacency" is a matrix that, for each (numeric) GO term id,
        // tells us which GO term are its "parents" following the required
        // relationships (is_a, part_of, ...) that were said to be parsed in the
        // program options
        // 2.- HSM computation       
        // 2.1.- Builds an HSM interfacer (to abstract the different HSMs)
        Object[] params = generateParameters(gti, matrixAxis);
        HSMInterfacer hsmi = buildsHSMInterfacer(params, new HashSet<GOTerm>(goIDsAsGOTerm), matrixAxis);
        hsmi.retrieveHSMinstance(this.hsmChoice, params);

        // 2.2.- Iterate and make the whole process for every desired ontology
        int loopVars[] = this.setLoopVars(dagChoice, logger);
        SolutionPrinter solutionPrinter = new SolutionPrinter(logger);
        for (int ontology = loopVars[0]; ontology < loopVars[1]; ontology++) {
            // for each ontology...
            // (a) compute HSM
            logger.showMessage("##### Computing HSM (" + new String[]{"BP", "MF", "CC"}[ontology] + ") #####");
            Matrix hsmResults;

            String genesRows[] = null;

            if (this.termWise) {
                // compute HSM term-wise
                hsmResults = hsmi.returnTermWiseResults(ontology);
            } else {
                hsmResults = hsmi.returnGeneWiseResults(ontology);
                if (this.geneIDs != null) {
                    genesRows = this.geneIDs;
                } else {
                    genesRows = hsmi.getComputedGenes();
                }
            }
            logger.log("HSM calculated");
            logger.showMemoryUsage();

            // (b) we print the results of the HSM to a file...            
            logger.showMessage("##### Printing HSM Results to File (" + new String[]{"BP", "MF", "CC"}[ontology] + ") #####");
            if (this.matrixStyle == ISM.MATRIX_STYLE || this.matrixStyle == ISM.BOTH_FILES) {
                solutionPrinter.printResultsToFile(ontology, hsmResults, matrixAxis, this.hsmFileName, this.notes, goIDsAsGOTerm, genesRows);
            }

            if (this.matrixStyle == ISM.TRIPLET_STYLE || this.matrixStyle == ISM.BOTH_FILES) {
                solutionPrinter.printeResultsToFileTripletStyle(ontology, hsmResults, matrixAxis, this.hsmFileName + "_triplet", this.notes, goIDsAsGOTerm, genesRows);
            }

            // (c) if we are to compute an ISM...
            if (this.isIsmToBeComputed) {
                hsmResults = hsmi.getOriginalCachedMatrix();
                // (d) we compute it

                ISMInterfacer ism = new ISMInterfacer();

                logger.showMessage("##### Computing ISM (" + new String[]{"BP", "MF", "CC"}[ontology] + ") #####");
                Matrix ismResults;
                if (this.termWise) {
                    // compute ISM term-wise
                    ismResults = ism.getISMs(matrixAxis, hsmResults, gti.getResults(), goIDsAsGOTerm, chosenRelations, dagChoice, ontology, logger);
                } else {
                    // compute ISM gene-wise
                    ismResults = ism.getGeneISMs(matrixAxis, hsmResults, gti.getResults(), chosenRelations, dagChoice, ontology, logger, this.weightedJaccard, this.geneIDs, hsmi.getComputedGenes());
                }

                // and we print the results of the HSM to a file...            
                logger.showMessage("##### Printing ISM Results to File (" + new String[]{"BP", "MF", "CC"}[ontology] + ") #####");
                logger.showMemoryUsage();
                if (this.matrixStyle == ISM.MATRIX_STYLE || this.matrixStyle == ISM.BOTH_FILES) {
                    solutionPrinter.printResultsToFile(ontology, ismResults, matrixAxis, this.ismFileName, this.notes, goIDsAsGOTerm, genesRows);
                }
                if (this.matrixStyle == ISM.TRIPLET_STYLE || this.matrixStyle == ISM.BOTH_FILES) {
                    solutionPrinter.printeResultsToFileTripletStyle(ontology, ismResults, matrixAxis, this.ismFileName + "_triplet", this.notes, goIDsAsGOTerm, genesRows);
                }
            }

        }
    }

    /**
     * Builds the HSM interfacer
     *
     * @param params set of parameters passed to the HSM interfacer
     * @return the HSMInterfacer object which has been built
     */
    private HSMInterfacer buildsHSMInterfacer(Object[] params, Set<GOTerm> targets, GOTerm[][] matrixAxis)
            throws IOException {
        //initialise the HSM interfacer

        HSMInterfacer hsmi = new HSMInterfacer(logger, targets, matrixAxis, this.geneIDs);
        logger.log("HSM_Interfacer instantiated");

        //initialise an instance of the required HSM method
        return hsmi;
    }

    private Object[] generateParameters(GOtreeInterfacer gti, GOTerm[][] matrixAxis) {
        //Create array of all GO terms
        GOTerm[] allterms = new GOTerm[matrixAxis[0].length + matrixAxis[1].length + matrixAxis[2].length];
        int counter = 0;
        for (int i = 0; i < 3; i++) {
            for (GOTerm go : matrixAxis[i]) {
                allterms[counter] = go;
                counter++;
            }
        }
        //set parameters for the HSMs
        Object[] params = new Object[]{allterms, gti.getGeneIDs(), gti.getGoIdsByGene(), matrixAxis, this.geneIDs, gti.getResults(), this.chosenRelations, this.logger}; //ALTERABLE PARAMETERS
        return params;
    }

    /**
     * Checks the choice of ontology; 'dagChoice' is acceptable, the 'log'
     * boolean variable refers to if a log file is being written & whether to
     * write any necessary information to it. This same use applies to the same
     * parameter in other methods of this class
     */
    private int[] setLoopVars(String dagChoice, TinyLogger logger) throws IOException {
        int[] loopVars = new int[2];
        if (dagChoice.compareToIgnoreCase("all") == 0) {
            loopVars[0] = 0;
            loopVars[1] = 3;
        } else if (dagChoice.compareToIgnoreCase("bp") == 0) {
            loopVars[0] = 0;
            loopVars[1] = 1;
        } else if (dagChoice.compareToIgnoreCase("mf") == 0) {
            loopVars[0] = 1;
            loopVars[1] = 2;
        } else if (dagChoice.compareToIgnoreCase("cc") == 0) {
            loopVars[0] = 2;
            loopVars[1] = 3;
        } else {
            logger.logAndCloseWriter("#######ERROR: Choice of Ontology Invalid");
            System.err.println("ERROR: Choice of Ontology Invalid");
            System.exit(-1);
        }
        return loopVars;
    }

    private void farewell() throws IOException {
        ui.farewell();
        this.logger.log("All printing done");
        this.logger.logAndCloseWriter("Finished");
    }
}
