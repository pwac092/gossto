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

import java.io.*;
import java.util.*;
import util.TinyLogger;

/**
 *
 * @author Samuel Heron
 */
//validates certain inputs, prints similarity matrices and creates the log file
public class IoValidation {

    private static TinyLogger logger;

    //Instantiates loop control variables
    public IoValidation() throws FileNotFoundException {
    }

    public static void setLogger(TinyLogger logga) {
        logger = logga;
    }

    //Checks the choice of ontology; 'dagChoice' is acceptable, the 'log' boolean variable refers to if a log file is being written & whether to write any necessary
    //information to it. This same use applies to the same parameter in other methods of this class 
    public static String validateDagChoice(String dagChoice) throws IOException {
        String choice = dagChoice.trim().toLowerCase();
        if (!Arrays.asList(new String[]{"all", "bp", "mf", "cc"}).contains(choice)) {
            logger.logAndCloseWriter("#######ERROR: Choice of Ontology Invalid");
            System.err.println("ERROR: Choice of Ontology Invalid");
            System.exit(-1);
        }
        return choice;
    }

    //Checks the choice of GO relations; 'relations', are acceptable
    public static void validateRelations(String[] relations) throws IOException {
        String acceptable_rels[] = {"is_a", "part_of", "regulates",
            "positively_regulates", "negatively_regulates", "has_part"};

        Set<String> acceptable = new HashSet<String>(Arrays.asList(acceptable_rels));
        for (String rel : relations) {
            if (!acceptable.contains(rel)) {
                logger.logAndCloseWriter("#########ERROR: Relation: " + rel + ", is not acceptable or does not exist");
                System.err.println("ERROR: Relation: " + rel + ", is not acceptable or does not exist");
                System.exit(-1);
            }
        }
    }

    //Checks the choice of evidence codes; 'codes', are acceptable
    public static void validateEvidenceCodes(String[] codes) throws IOException {

        String acceptable_codes[] = {"EXP", "IDA", "IPI", "IMP", "IGI",
            "IEP", "TAS", "IC", "ISS", "ISO", "ISA", "ISM", "IGC",
            "IBA", "IBD", "IKR", "IRD", "RCA", "NAS", "ND", "IEA", "ALL"};

        Set<String> acceptable = new HashSet<String>(Arrays.asList(acceptable_codes));

        for (String code : codes) {
            if (!acceptable.contains(code.toUpperCase())) {
                logger.logAndCloseWriter("#######ERROR: Evidence code: " + code + ", is not acceptable or does not exist.");
                System.err.println("ERROR: Evidence code: " + code + ", is not acceptable or does not exist.");
                System.exit(-1);
            }
        }
    }

    //Checks that supplied file paths; 'oboFile' & 'goaFile', exist
    public static void validateFilePaths(String oboFile, String goaFile) throws IOException {
        if (!(new File(oboFile).isFile())) {
            logger.logAndCloseWriter("ERROR: No OBO file at specified filepath");
            System.err.println("########ERROR: No OBO file at specified filepath");
            System.exit(-1);
        }
        if (!(new File(goaFile).isFile())) {
            logger.logAndCloseWriter("##########ERROR: No GOA file at specified filepath");
            System.err.println("ERROR: No GOA file at specified filepath");
            System.exit(-1);
        }
    }

    //Checks that supplied output location; 'location', exists
    public static void validateOutputLocation(String location) throws IOException {

        if (location.contains("/")) {
            for (int i = location.length() - 1; i > 1; i--) //finding the last '/', -1 as length counts 0 index
            {
                if (location.charAt(i) == '/') {
                    try {
                        File output = new File(location.substring(0, i));
                        if (output.isDirectory() && output.canWrite() && output.canRead() && output.isAbsolute()) {
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
}
