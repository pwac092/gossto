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
import Jama.Matrix;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.TinyLogger;

/**
 *
 * @author aeromero
 */
public class SolutionPrinter {

    TinyLogger logger;
    private static final String ontologies[] = {"Biological process", "Molecular function", "Cellular Component"};

    public SolutionPrinter(TinyLogger logger) {
        this.logger = logger;
    }

    /**
     * @param ontology identifier of the ontology (0, 1 or 2)
     * @param matrix matrix with the computed results (already 'reduced' matrix)
     * @param axis
     * @param outputName name of the output file
     * @param notes
     * @param geneIDs
     * @param targetGoIDs
     * @throws java.io.IOException
     */
    public void printResultsToFile(int ontology, Matrix matrix, GOTerm[][] axis, String outputName, ArrayList<String> notes, ArrayList<GOTerm> targetGoIDs, String[] geneIDs) throws IOException {
        //re-validate file path:
        IoValidation.validateOutputLocation(outputName);
        File outputFileName = getOutputFileName(ontology, outputName);

        /* Filter null or very small matrices */
        if (matrix == null) {
            printMessageNotEnoughAnnotations(outputFileName);
            return;
        }
        //printing the matrices value by value
        final int n = matrix.getRowDimension();
        final int m = matrix.getColumnDimension();
        if (n == 1) {
            logger.showMessage("  ERROR: the specified user restrictions leaves a 1x1 matrix,");
            logger.showMessage("  which will not be printed.");
            return;
        }

        String[] rowIdentifiers = getRowIdentifiers(n, targetGoIDs, getGOIds(targetGoIDs, ontology), geneIDs, axis, ontology);

        try {
            logger.showMessage("Printing results for Ontology : " + ontologies[ontology]);

            BufferedWriter out = new BufferedWriter(new FileWriter(outputFileName), 32768);
            // write notes
            for (String note : notes) {
                out.write("! ");
                out.write(note);
                out.newLine();
            }

            logger.showMessage("Notes printed");
            // write first line
            for (String row : rowIdentifiers) {
                out.write(row);
                out.write("\t");
            }

            out.newLine();
            for (int i = 0; i < n; i++) {
                out.write(rowIdentifiers[i]);
                for (int j = 0; j < m; j++) {
                    //Checks the size of the similarity value, if smaller than 0.001 or greater than -0.001 then some validation needs to take place to print them properly 
                    //double checkVal = matrix.getEntry(i, j);
                    //String outputFormat = this.myFormat(checkVal);
                    out.write("\t");
                    out.write("" + matrix.get(i, j));
                }

                out.newLine();
            }
            out.close();
            logger.log("Printing complete; Output File: " + outputFileName);
            System.out.println("Printing COMPLETE; Output File: " + outputFileName);

        } catch (java.lang.OutOfMemoryError oome) {
            logger.logAndCloseWriter("############## ERROR: Out of memory Error of type: " + oome.getMessage());
            System.err.println("ERROR: Java has run out of memory. Memory Type: " + oome.getMessage());
            System.exit(-1);
        }
    }

    private Set<Integer> getGOIds(ArrayList<GOTerm> targetGoIDs, int ontology) {
        Map<Integer, Set<Integer>> goTermIds = new HashMap<Integer, Set<Integer>>();
        goTermIds.put(0, getTermsForOntology(targetGoIDs, "biological_process"));
        goTermIds.put(1, getTermsForOntology(targetGoIDs, "molecular_function"));
        goTermIds.put(2, getTermsForOntology(targetGoIDs, "cellular_component"));
        Set<Integer> goIds = goTermIds.get(ontology);
        return goIds;
    }

    public void printeResultsToFileTripletStyle(int ontology, Matrix matrix, GOTerm[][] axis, String outputName, ArrayList<String> notes, ArrayList<GOTerm> targetGoIDs, String[] geneIDs) throws IOException {
        //re-validate file path:
        IoValidation.validateOutputLocation(outputName);
        //Creates output file(s) where specified
        File outputFileName = getOutputFileName(ontology, outputName);
        if (matrix == null) {
            printMessageNotEnoughAnnotations(outputFileName);
            return;
        }
        //printing the matrices value by value
        final int n = matrix.getRowDimension();
        final int m = matrix.getColumnDimension();
        if (n == 1) {
            logger.showMessage("  ERROR: the specified user restrictions leaves a 1x1 matrix,");
            logger.showMessage("  which will not be printed.");
            return;
        }

        try {
            logger.showMessage("Printing results for Ontology : " + ontologies[ontology]);

            String[] rowIdentifiers = getRowIdentifiers(n, targetGoIDs, getGOIds(targetGoIDs, ontology), geneIDs, axis, ontology);

            logger.showMessage("Printing contents: " + n + " " + m);
            BufferedWriter out = new BufferedWriter(new FileWriter(outputFileName), 32768);
            for (int i = 0; i < n - 1; i++) {
                for (int j = i + 1; j < n; j++) {
                    out.write(rowIdentifiers[i]);
                    out.write("\t");
                    out.write(rowIdentifiers[j]);
                    out.write("\t" + matrix.get(i, j));
                    out.newLine();
                }
            }
            out.close();
            logger.log("Printing complete; Output File: " + outputFileName);
            System.out.println("Printing COMPLETE; Output File: " + outputFileName);
        } catch (java.lang.OutOfMemoryError oome) {
            logger.logAndCloseWriter("############## ERROR: Out of memory Error of type: " + oome.getMessage());
            System.err.println("ERROR: Java has run out of memory. Memory Type: " + oome.getMessage());
            System.exit(-1);
        }
    }

    private String[] getRowIdentifiers(final int n, ArrayList<GOTerm> targetGoIDs, Set<Integer> goIds, String[] geneIDs, GOTerm[][] axis, int ontology) {
        String[] rowIdentifiers = new String[n];
        int ind = 0;
        if (targetGoIDs != null && !targetGoIDs.isEmpty()) {
            for (GOTerm term : targetGoIDs) {
                if (goIds.contains(term.getNumericId())) {
                    rowIdentifiers[ind++] = term.getGOid();
                }
            }
        } else if (geneIDs != null) {
            System.arraycopy(geneIDs, 0, rowIdentifiers, 0, n);
        } else {
            for (GOTerm term : axis[ontology]) {
                rowIdentifiers[ind++] = term.getGOid();
            }
        }
        return rowIdentifiers;
    }

    private void printMessageNotEnoughAnnotations(File outputFileName) throws IOException {
        logger.showMessage("\n  ERROR: the specified annotation file has not enough data to produce a semantic similarity matrix for this ontology.\n");
        logger.showMessage("  This may be due to one of the following facts:\n\n");
        logger.showMessage("  - There are not enough annotations (possibly none) with the specified evidence codes (try adding more, especially IEA).");
        logger.showMessage("  - The organism has annotations, but the subset of GO terms you selected contains no annotation (try being less restrictive).");
        logger.showMessage("  - There are not enough genes (possibly none) annotated to valid entries (try selecting a different organism).");
        logger.showMessage("  - The organism has annotations, but the subset of genes you selected contains no associated GO terms (try being less restrictive).");
        logger.showMessage("  We cannot do much for solving your problem. We recommend you to select a better organism,");
        logger.showMessage("  or maybe download GOssTo in its Java version from our webpage ( http://www.paccanarolab.org/gossto/ ),");
        logger.showMessage("  and tweak it a bit to cope with this annotation file.\n\n");

        BufferedWriter out = new BufferedWriter(new FileWriter(outputFileName));
        out.write("ERROR: the specified annotation file has not enough data to produce a semantic similarity matrix for this ontology.");
        out.newLine();
        out.newLine();
        out.write("This may be due to one of the following facts:");
        out.newLine();
        out.write("- There are not enough annotations (possibly none) with the specified evidence codes (try adding more, especially IEA).");
        out.newLine();
        out.write("- The organism has annotations, but the subset of GO terms you selected contains no annotation (try being less restrictive).");
        out.newLine();
        out.write("- There are not enough genes (possibly none) annotated to valid entries (try selecting a different organism).");
        out.newLine();
        out.write("- The organism has annotations, but the subset of genes you selected contains no associated GO terms (try being less restrictive).");
        out.newLine();
        out.write("We cannot do much for solving your problem. We recommend you to select a better organism,");
        out.newLine();
        out.write("or maybe download GOssTo in its Java version from our webpage ( http://www.paccanarolab.org/gossto/ ),");
        out.write(" and tweak it a bit to cope with this annotation file.");
        out.close();
    }

    private File getOutputFileName(int ontology, String outputName) {
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

        String ext[] = {"_BP.txt", "_MF.txt", "_CC.txt"};

        if (addExtension) {
            return new File(outputName + ext[ontology]);
        } else {
            return new File(outputName);
        }

    }

    private Set<Integer> getTermsForOntology(List<GOTerm> terms, String ontologyName) {
        HashSet<Integer> myOntologies = new HashSet<Integer>();
        if (terms != null) {
            for (GOTerm term : terms) {
                if (term.getOntology().getName().equals(ontologyName)) {
                    myOntologies.add(term.getNumericId());
                }
            }
        }
        return myOntologies;
    }

}
