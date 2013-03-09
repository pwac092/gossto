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

import GOtree.GOTerm;
import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.apache.commons.math3.linear.RealMatrix;
import util.TinyLogger;

/**
 *
 * @author Samuel Heron
 */
//validates certain inputs, prints similarity matrices and creates the log file
public class IoValidation {
    
    private static final String ontologies[] = {"Biological process", "Molecular function", "Cellular Component"};
    private int[] loopVars; //Used to control calculations/printing based upon the choice of ontology to save time
    private TinyLogger logger;

    //Instantiates loop control variables
    public IoValidation() throws FileNotFoundException {
        loopVars = new int[2];
    }
    
    public void setLogger(TinyLogger logger) {
        this.logger = logger;
    }

    //returns the loop control variables
    public int[] getLoopVariables() {
        return this.loopVars;
    }

    //Checks the choice of ontology; 'dagChoice' is acceptable, the 'log' boolean variable refers to if a log file is being written & whether to write any necessary
    //information to it. This same use applies to the same parameter in other methods of this class 
    public String validateDagChoice(String dagChoice) throws IOException {
        if (dagChoice.compareToIgnoreCase("all") == 0) {
            dagChoice = "all";
            loopVars[0] = 0;
            loopVars[1] = 3;
        } else if (dagChoice.compareToIgnoreCase("bp") == 0) {
            dagChoice = "bp";
            loopVars[0] = 0;
            loopVars[1] = 1;
        } else if (dagChoice.compareToIgnoreCase("mf") == 0) {
            dagChoice = "mf";
            loopVars[0] = 1;
            loopVars[1] = 2;
        } else if (dagChoice.compareToIgnoreCase("cc") == 0) {
            dagChoice = "cc";
            loopVars[0] = 2;
            loopVars[1] = 3;
        } else {
            logger.logAndCloseWriter("#######ERROR: Choice of Ontology Invalid");
            System.err.println("ERROR: Choice of Ontology Invalid");
            System.exit(-1);
        }
        return dagChoice;
    }

    //Checks the choice of GO relations; 'relations', are acceptable
    public void validateRelations(String[] relations) throws IOException {
        Set<String> acceptable = new HashSet<String>();
        acceptable.add("is_a");
        acceptable.add("part_of");
        acceptable.add("regulates");
        acceptable.add("positively_regulates");
        acceptable.add("negatively_regulates");
        acceptable.add("has_part");
        for (String rel : relations) {
            if (acceptable.contains(rel) == false) {
                logger.logAndCloseWriter("#########ERROR: Relation: " + rel + ", is not acceptable or does not exist");
                System.err.println("ERROR: Relation: " + rel + ", is not acceptable or does not exist");
                System.exit(-1);
            }
        }
    }

    //Checks the choice of evidence codes; 'codes', are acceptable
    public void validateEvidenceCodes(String[] codes) throws IOException {
        Set<String> acceptable = new HashSet<String>();
        acceptable.add("EXP");
        acceptable.add("IDA");
        acceptable.add("IPI");
        acceptable.add("IMP");
        acceptable.add("IGI");
        acceptable.add("IEP");
        acceptable.add("TAS");
        acceptable.add("IC");
        acceptable.add("ISS");
        acceptable.add("ISO");
        acceptable.add("ISA");
        acceptable.add("ISM");
        acceptable.add("IGC");
        acceptable.add("IBA");
        acceptable.add("IBD");
        acceptable.add("IKR");
        acceptable.add("IRD");
        acceptable.add("RCA");
        acceptable.add("NAS");
        acceptable.add("ND");
        acceptable.add("IEA");
        acceptable.add("ALL");
        for (String code : codes) {
            if (acceptable.contains(code.toUpperCase()) == false) {
                logger.logAndCloseWriter("#######ERROR: Evidence code: " + code + ", is not acceptable or does not exist.");
                System.err.println("ERROR: Evidence code: " + code + ", is not acceptable or does not exist.");
                System.exit(-1);
            }
        }
    }

    //Checks that supplied file paths; 'oboFile' & 'goaFile', exist
    public void validateFilePaths(String oboFile, String goaFile) throws IOException {
        if (new File(oboFile).isFile() == false) {
            logger.logAndCloseWriter("ERROR: No OBO file at specified filepath");
            System.err.println("########ERROR: No OBO file at specified filepath");
            System.exit(-1);
        }
        if (new File(goaFile).isFile() == false) {
            logger.logAndCloseWriter("##########ERROR: No GOA file at specified filepath");
            System.err.println("ERROR: No GOA file at specified filepath");
            System.exit(-1);
        }
    }

    //Checks that supplied output location; 'location', exists
    public void validateOutputLocation(String location) throws IOException {
        
        if (location.contains("/") == true) {
            for (int i = location.length() - 1; i > 1; i--) //finding the last '/', -1 as length counts 0 index
            {
                if (location.charAt(i) == '/') {
                    try {
                        File output = new File(location.substring(0, i));
                        if (output.isDirectory() == true && output.canWrite() == true && output.canRead() == true && output.isAbsolute() == true) {
                            return;
                        } else {
                            logger.logAndCloseWriter("########ERROR: Invalid output location");
                            System.err.println("ERROR: Invalid output location");
                            System.exit(-1);
                        }
                    } catch (Exception e) {
                        logger.logAndCloseWriter("########ERROR: Invalid output location: " + e.getMessage());
                        System.err.println("ERROR: Invalid output location: " + e.getMessage());
                        System.exit(-1);
                    }
                }
            }
        }
    }
    
    private Set<Integer> getTermsForOntology(List<GOTerm> terms, String ontologyName) {
        HashSet<Integer> myOntologies = new HashSet<Integer>();
        if (null != terms) {
            for (GOTerm term : terms) {
                if (term.getOntology().getName().equals(ontologyName)) {
                    myOntologies.add(term.getNumericId());
                }
            }
        }
        return myOntologies;
    }
    
    public void printResultsToFile(int ontology, RealMatrix matrix, GOTerm[][] axis, String outputName, ArrayList<String> notes, ArrayList<GOTerm> targetGoIDs, String[] geneIDs) throws IOException {
        //re-validate file path:
        validateOutputLocation(outputName);
        //Creates output file(s) where specified
        File outputFileName;

        //check the output names to make sure we print a friendly name.
        //0. we are looking for .something. this something, would ussually
        //be three characters long. But, since POSIX permits any length of names.
        //however, since people might choose to use dots in the middle of the 
        //file name (weird, but who are we to judge..) we will just accept those
        //names and any and all names with a dot anywhere, we'll consider 
        //as full names, and will no append anything to them. Otherwise, 
        //just to make things easier, we will append .txt termination. Why txt?
        //Well, for someone not that experienced, this is easily recongnisable
        //and for someone with a bit more experience, a double click will
        //not mean opening some spreadsheet program. Long comment for a simple
        //feature, but..
        
        Pattern p = Pattern.compile("\\.*");
        Matcher matcher = p.matcher(outputName);
        boolean addExtension = !matcher.matches();
        
        
        switch (ontology) {
            case 0:
                if (addExtension) {
                    outputFileName = new File(outputName + "_BP.txt");
                } else {
                    outputFileName = new File(outputName);
                }
                break;
            case 1:
                if (addExtension) {
                    outputFileName = new File(outputName + "_MF.txt");
                } else {
                    outputFileName = new File(outputName);
                }
                break;
            default:
                if (addExtension) {
                    outputFileName = new File(outputName + "_CC.txt");
                } else {
                    outputFileName = new File(outputName);
                }
        }
        
        Map<Integer, Set<Integer>> goTermIds = new HashMap<Integer, Set<Integer>>();
        goTermIds.put(0, getTermsForOntology(targetGoIDs, "biological_process"));
        goTermIds.put(1, getTermsForOntology(targetGoIDs, "molecular_function"));
        goTermIds.put(2, getTermsForOntology(targetGoIDs, "cellular_component"));
        
        Set<Integer> goIds = goTermIds.get(ontology);
        try {
            if (matrix != null) {
                logger.showMessage("Printing results for Ontology : " + ontologies[ontology]);
                
                BufferedWriter out = new BufferedWriter(new FileWriter(outputFileName), 32768);
                //###file notes###
                for (String note : notes) {
                    out.write("! ");
                    out.write(note);
                    out.newLine();
                }
                
                logger.showMessage("Notes printed");
                //put together the top axis
                String[] temp = new String[matrix.getRowDimension()];
                logger.showMessage("Row dimension of the matrix: " + matrix.getRowDimension());
                int ind = 0;
                if (targetGoIDs != null && !targetGoIDs.isEmpty()) {
                    int termsWritten = 0;
                    for (GOTerm term : targetGoIDs) {
                        if (goIds.contains(term.getNumericId())) {
                            if (termsWritten > 0) {
                                out.write("\t");
                            }
                            temp[ind] = term.getGOid();
                            out.write(temp[ind]);
                            ind++;
                            termsWritten++;
                            //lineToWrite += term.getGOid() + ",";
                        }
                    }
                } else if (geneIDs != null) {
                    int genesWritten = 0;
                    
                    for (String term : geneIDs) {
                        if (genesWritten > 0) {
                            out.write("\t");
                        }
                        out.write(term);
                        temp[ind] = term;
                        ind++;
                        ++genesWritten;
                    }
                } else {
                    int termsWritten = 0;
                    
                    for (GOTerm term : axis[ontology]) {
                        if (termsWritten > 0) {
                            out.write("\t");
                        }
                        out.write(term.getGOid() + "");
                        temp[ind] = term.getGOid();
                        ind++;
                        
                        ++termsWritten;
                    }
                }
                //out.write(lineToWrite.substring(0, lineToWrite.length() - 1));
                out.newLine();
                //lineToWrite = "";

                //printing the matrices value by value
                final int n = matrix.getRowDimension();
                final int m = matrix.getColumnDimension();
                
                logger.showMessage("Printing contents: " + n + " " + m);
                
                for (int i = 0; i < n; i++) {
                    out.write(temp[i]);
                    out.write("\t");
                    
                    for (int j = 0; j < m; j++) {
                        //Checks the size of the similarity value, if smaller than 0.001 or greater than -0.001 then some validation needs to take place to print them properly 
                        //double checkVal = matrix.getEntry(i, j);
                        //String outputFormat = this.myFormat(checkVal);

                        out.write("" + (float) matrix.getEntry(i, j));
                        if (j != m - 1) {
                            out.write("\t");
                        }
                    }
                    
                    out.newLine();
                }
                out.close();
                logger.log("Printing complete; Output File: " + outputFileName);
                System.out.println("Printing COMPLETE; Output File: " + outputFileName);
            } else {
                
                logger.showMessage("\n  ERROR: the specified annotation file has not enough data to produce a semantic similarity matrix for this ontology.\n");
                logger.showMessage("  This may be due to one of the following facts:\n\n");
                logger.showMessage("  - There are not many annotations (possibly none) with the specified evidence codes (try adding more, especially IEA).");
                logger.showMessage("  - The organism has annotations, but the subset of GO terms you selected contains no annotation (try being less restrictive).");
                logger.showMessage("  - There are not many genes (possibly none) annotated to valid entries (try selecting a different organism).");
                logger.showMessage("  - The organism has annotations, but the subset of genes you selected contains no associated GO terms (try being less restrictive).");
                logger.showMessage("  We cannot do much for solving your problem. We recommend you to select a better organism,");
                logger.showMessage("  or maybe download GOssTo in its Java version from our webpage ( http://www.paccanarolab.org/gossto/ ),");
                logger.showMessage("  and tweak it a bit to cope with this annotation file.\n\n");
                
                
                   BufferedWriter out = new BufferedWriter(new FileWriter(outputFileName), 32768);                
                out.write("ERROR: the specified annotation file has not enough data to produce a semantic similarity matrix for this ontology.");
                out.newLine();
                out.newLine();
                out.write("This may be due to one of the following facts:");
                out.newLine();
                out.write("- There are not many annotations (possibly none) with the specified evidence codes (try adding more, especially IEA).");
                out.newLine();
                out.write("- The organism has annotations, but the subset of GO terms you selected contains no annotation (try being less restrictive).");
                out.newLine();                
                out.write("- There are not many genes (possibly none) annotated to valid entries (try selecting a different organism).");
                out.newLine();
                out.write("- The organism has annotations, but the subset of genes you selected contains no associated GO terms (try being less restrictive).");
                out.newLine();
                out.write("We cannot do much for solving your problem. We recommend you to select a better organism,");
                out.newLine();
                out.write("or maybe download GOssTo in its Java version from our webpage ( http://www.paccanarolab.org/gossto/ ),");
                out.write(" and tweak it a bit to cope with this annotation file.");
                out.close();
                
            }
        } catch (java.lang.OutOfMemoryError oome) {
            logger.logAndCloseWriter("############## ERROR: Out of memory Error of type: " + oome.getMessage());
            System.err.println("ERROR: Java has run out of memory. Memory Type: " + oome.getMessage());
            System.exit(-1);
        }
    }
}
