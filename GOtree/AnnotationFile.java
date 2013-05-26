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
package GOtree;

/*
 * To change this template, choose Tools | Templates and open the template  n
 * the editor.
 */
import AnnotationFileUtilities.AnnotationFileStrategy;
import AnnotationFileUtilities.AnnotationFileWithExternalProteinIds;
import AnnotationFileUtilities.AnnotationFileWithoutExternalProteinIds;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/**
 *
 * @author aeromero
 */
public class AnnotationFile {

    public static final int USE_ALL = 0;
    public static final int USE_ALL_BUT_IEA = 1;
    public static final int USE_EXP = 2;
    /**
     * Experimental codes used in the Gene Ontology Annotation
     */
    private final static Set<String> experimentalCodes = new HashSet<String>(Arrays.asList(new String[]{"EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "TAS", "IC"}));
    /**
     * All evidence codes used in the Gene Ontology Annotation, except IEA
     */
    private final static Set<String> allCodesButIEA = new HashSet<String>(Arrays.asList(new String[]{"EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "TAS", "IC", "ISS", "ISO", "ISA", "ISM", "IGC", "IBA", "IBD", "IKR", "IRD", "RCA", "NAS", "ND"}));
    /**
     * All evidence codes used in the Gene Ontology Annotation
     */
    private final static Set<String> allCodes = new HashSet<String>(Arrays.asList(new String[]{"EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "TAS", "IC", "ISS", "ISO", "ISA", "ISM", "IGC", "IBA", "IBD", "IKR", "IRD", "RCA", "NAS", "ND", "IEA", "NR"}));
    
    /**
     * Identifiers of the proteome which match this annotation file
     */
    private Set<String> proteinIdentifiers;
    /**
     * Strategy we are using to process this annotation file, whether we are
     * processing all protein identifiers, or we are just considering a set of
     * identifiers specified by hand
     */
    AnnotationFileStrategy strategy;

    public AnnotationFile() {
        this.proteinIdentifiers = new HashSet<String>();
    }

    public void setProteinIdentifiers(Set<String> proteinIds) {
        this.proteinIdentifiers = new HashSet<String>();
        this.proteinIdentifiers.addAll(proteinIds);
    }

    private Set<String> getSetOfEvidenceCodes(int evidenceCodes) {
        assert (evidenceCodes == AnnotationFile.USE_ALL || evidenceCodes == AnnotationFile.USE_ALL_BUT_IEA || evidenceCodes == AnnotationFile.USE_EXP);

        Set<String> setOfEvidenceCodes = null;
        switch (evidenceCodes) {
            case AnnotationFile.USE_ALL:
                setOfEvidenceCodes = AnnotationFile.allCodes;
                break;
            case AnnotationFile.USE_ALL_BUT_IEA:
                setOfEvidenceCodes = AnnotationFile.allCodesButIEA;
                break;
            case AnnotationFile.USE_EXP:
                setOfEvidenceCodes = AnnotationFile.experimentalCodes;
                break;
        }
        return setOfEvidenceCodes;
    }

    private Assignment readFromFile(String fileName, final Set<String> validEvicenceCodes) throws FileNotFoundException, IOException {
        BufferedReader reader = new BufferedReader(new FileReader(fileName));
        Assignment output = new Assignment();

        while (reader.ready()) {
            String line = reader.readLine();

            if (line.startsWith("!")) {
                continue;
            }

            final String fields[] = line.split("\\t");
            final String evidenceCode = fields[6];

            if (!validEvicenceCodes.contains(evidenceCode)) {
                continue;
            }

            if (this.strategy.proteinIdIsValid(fields)) {
                final String proteinId = this.strategy.getProteinId(fields);
                final String goTerm = fields[4];
                output.setValue(proteinId, goTerm, 1.0);
            }
        }
        reader.close();
        return output;
    }

    public Set<String> getUnparsedProteins() {
        return this.strategy.getSetOfUnParsedProteins();
    }

    public static void useUniProtIds(boolean val) {
        AnnotationFileWithoutExternalProteinIds.setUseUniProtIds(val);
    }

    public Assignment readAnnotationFileWithProteinIdentifiers(String fileName, Set<String> proteinIds, int evidenceCodes) throws IOException, FileNotFoundException {

        Set<String> setOfEvidenceCodes = this.getSetOfEvidenceCodes(evidenceCodes);
        this.strategy = new AnnotationFileWithExternalProteinIds();
        this.strategy.setProteinIdentifiers(proteinIds);
        return this.readFromFile(fileName, setOfEvidenceCodes);
    }

    // reads annotation file with all evidence codes
    public Assignment readAnnotationFile(String fileName) throws IOException, FileNotFoundException {
        this.strategy = new AnnotationFileWithoutExternalProteinIds();
        return this.readFromFile(fileName, allCodes);
    }

    //I edited this so that the evidence codes are supplied directly rather than an integer choice relayed to the getSetOfEvidenceCodes method
    public Assignment readAnnotationFile(String fileName, String[] evidenceCodes) throws IOException, FileNotFoundException {
        Set<String> setOfEvidenceCodes = new HashSet<String>();
        setOfEvidenceCodes.addAll(Arrays.asList(evidenceCodes));
        this.strategy = new AnnotationFileWithoutExternalProteinIds();
        return this.readFromFile(fileName, setOfEvidenceCodes);
    }

    /*public static void main(String[] args) throws FileNotFoundException, IOException {
    AnnotationFile file = new AnnotationFile();
    
    Assignment ass = null;
    
    String names[] = {"gene_association.PAMGO_Oomycetes"};//{"goa1", "goa2", "goa3"};
    
    for (String name : names) {
    System.err.println("Reading file " + name + " (all. codes)");
    ass = file.readAnnotationFile(name, AnnotationFile.USE_ALL);
    System.err.println("=====================================================================");
    }
    
    for (String name : names) {
    System.err.println("Reading file " + name + " (exp. codes)");
    ass = file.readAnnotationFile(name, AnnotationFile.USE_EXP);
    }
    }*/
}
