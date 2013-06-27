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

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import util.TinyLogger;

/**
 *
 * @author aeromero
 */
public abstract class ParameterValidator {

    protected String[] chosenRelations;
    protected String[] evidenceCodes;
    protected String dagChoice;
    protected String ismFileName;
    protected String hsmFileName;
    protected String oboFile;
    protected boolean ismChoice;
    protected boolean termWise;
    protected String[] goIDs;
    protected String[] geneIDs;
    protected String hsmChoice;
    protected String goaFile;
    protected ArrayList<String> notes;
    protected boolean ownHSM;
    protected boolean weightedJaccard;
    protected int matrixStyle;
    protected boolean useUniProtIds;

    protected ParameterValidator() {
        oboFile = "";
        goaFile = "";
        hsmChoice = "";
        dagChoice = "";
        hsmFileName = "";
        ismFileName = "";

        chosenRelations = null;
        geneIDs = null;
        goIDs = null;
        evidenceCodes = null;
        //Stores the execution parameters to write them at the top of the output files
        notes = new ArrayList<String>();
        //boolean values relating to input values and choices
        ismChoice = true;
        termWise = true;
        weightedJaccard = false;
        ownHSM = false;
        this.useUniProtIds = true;
        this.matrixStyle = ISM.MATRIX_STYLE;
    }

    public abstract void validate(IoValidation validate, TinyLogger logger) throws FileNotFoundException, IOException;

    public String[] getChosenRelations() {
        return chosenRelations;
    }

    public String[] getEvidenceCodes() {
        return evidenceCodes;
    }

    public String getDagChoice() {
        return dagChoice;
    }

    public String getIsmFileName() {
        return ismFileName;
    }

    public String getHsmFileName() {
        return hsmFileName;
    }

    public String getOboFile() {
        return oboFile;
    }

    public boolean isIsmChoice() {
        return ismChoice;
    }

    public boolean isTermWise() {
        return termWise;
    }

    public String[] getGoIDs() {
        return goIDs;
    }

    public String[] getGeneIDs() {
        return geneIDs;
    }

    public String getHsmChoice() {
        return hsmChoice;
    }

    public String getGoaFile() {
        return goaFile;
    }

    public ArrayList<String> getNotes() {
        return notes;
    }

    public boolean isOwnHSM() {
        return ownHSM;
    }

    public boolean isWeightedJaccard() {
        return weightedJaccard;
    }

    public int getMatrixStyle() {
        return matrixStyle;
    }

    public boolean isUseUniProtIds() {
        return useUniProtIds;
    }
    
    
}
